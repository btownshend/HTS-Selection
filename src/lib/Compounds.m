% Data structure to hold information about compounds located in mass spec runs
classdef Compounds < handle
  properties
    names;   % names{i} - Name of compound i (e.g. 'CDIV0051-B07')
    shortnames; % shortnames{i} - Short name of compound i (e.g. '51B07')
    sdf;     % sdf{i} - SDF data for compound i
    mztarget; % mztarget(i) - target m/z for compound i 
    files;   % files{j} - Mass spec filename j
    moles;    % moles(j) - Moles of each compound loaded in run j
    group;   % group{j} - name of group this file belongs to
    contains;  % contains(i,j) is true if we expected to find compound i in file j
             % For compounds uniquely located in a particular run:
    mz;    % mz(i,j) contains the observed m/z peak on compound i, from file j
    time;    % time(i,j) contains the elution time on compound i, from file j
    ic;    % ic(i,j) contains the total ion count for compound i, from file j
    nisomers;   % nisomers(i,j) is the number of other compounds that are within MZFUZZ in this file
    nunique;   % nunique(i,j) is the number of other compounds that are within MZFUZZ in this file and with unknown elution time
    numhits;   % numhits(i,j) is the number of possible hits (need to have exactly 1 for the hit to be used)
    multihits;
  end
  
  properties(Constant)
    MZFUZZ=0.003;
    TIMEFUZZ=50;   % in seconds
  end
  
  methods
    function obj=Compounds()
      obj.multihits={};
      obj.contains=false(0,0);
    end

    function ind=find(obj,name)
    % Find given name
    % Try both format stored in names (e.g. 'CDIV0051-B07') and short format (e.g. '51B07')
      if name(end-1)>='A' && name(end-2)<='H'
        % Add leading 0 to column
        name=[name(1:end-1),'0',name(end)];
      end
      ind=find(strcmp(name,obj.names));
      if length(ind)~=1
        ind=find(strcmp(name,obj.shortnames));
        if length(ind)~=1
          error('Compound %s not found',name);
        end
      end
    end
      
    function nindex=lookupName(obj,name, mztarget, sdf)
    % Find the index of compound by name, create if missing
      if ismember(name,obj.names)
        nindex=find(strcmp(name,obj.names));
      else
        obj.names{end+1,1}=name;
        nindex=length(obj.names);
        if nargin>=4
          obj.sdf{nindex,1}=sdf;
        else
          obj.sdf{nindex,1}=[];
        end
        if nargin>=3
          obj.mztarget(nindex,1)=mztarget;
        else
          obj.mztarget(nindex,1)=nan;
        end
        obj.mz(nindex,:)=nan;
        obj.time(nindex,:)=nan;
        obj.ic(nindex,:)=nan;
        obj.contains(nindex,:)=false;
        
        % Set short name
        [p,r,c]=obj.getposition(length(obj.names));
        obj.shortnames{length(obj.names)} = sprintf('%d%c%02d',str2num(p(5:end)),r+'A'-1,c);
      end
      if length(obj.shortnames) < length(obj.names)
        obj.shortnames=cell(size(obj.names));
        for i=1:length(obj.names)
          [p,r,c]=obj.getposition(i);
          obj.shortnames{i} = sprintf('%d%c%02d',str2num(p(5:end)),r+'A'-1,c);
        end
      end
    end

    function findex=lookupMS(obj,ms)
    % Find the index of a file by name, create if missing
      if ismember(ms.name,obj.files)
        findex=find(strcmp(ms.name,obj.files));
      else
        obj.files{1,end+1}=ms.name;
        findex=length(obj.files);
        obj.mz(:,findex)=nan;
        obj.time(:,findex)=nan;
        obj.ic(:,findex)=nan;
        obj.contains(:,findex)=false;
        obj.moles(findex)=ms.moles;
      end
    end

    function [p,r,c]=getposition(obj,i)
    % Get plate (str),row (int),col(int) of compound i
      sdf=obj.sdf{i};
      for ii=1:length(i)
        p=sdf.BATCH_PLATE;
        r=sdf.BATCH_WELL(1)-'A'+1;
        c=sscanf(sdf.BATCH_WELL(2:end),'%d');
      end
    end
    
    function id=checkComposition(obj,ms,varargin)
    % Check composition in mass spec file using current set of compounds
      defaults=struct('debug',false);  
      args=processargs(defaults,varargin);
      id=[];
      for i=1:length(obj.mztarget)
        meantime(i)=nanmean(obj.time(i,obj.contains(i,:)));
      end
      for i=1:length(obj.mztarget)
        mztarget=obj.mztarget(i);
        id=[id,struct('mztarget',mztarget,'desc','','findargs','','mz',nan,'time',nan,'ic',nan)];
        if isfinite(meantime(i))
          % Only look for compounds with a known elution time
          isomers=find(abs(meantime-meantime(i))<obj.TIMEFUZZ & abs(mztarget-obj.mztarget')<obj.MZFUZZ);
          id(i)=ms.findcompound(mztarget,'sdf',obj.sdf{i},'elutetime',meantime(i),'timetol',obj.TIMEFUZZ,'mztol',obj.MZFUZZ,'debug',args.debug);
          if isempty(id(i).ic)
            id(i).ic=0;
            id(i).time=nan;
            id(i).mz=nan;
          end
          if length(isomers)>1 && sum(id(i).ic)>0
            %fprintf('Isomer at compounds %s with ion count = %.0f\n', sprintf('%d,',isomers),sum(id(i).ic));
            id(i).ic=nan;   % Can't use it
          end
        end
      end
      
      for i=1:length(obj.mztarget)
        id(i).name=obj.names{i};
        id(i).relic=sum(id(i).ic)/nanmax(obj.ic(i,:));
      end
    end
    
    function [matic,id,refid]=plotComposition(obj,ms,varargin)
      defaults=struct('debug',false,'thresh',0.02,'ref',[]);  
      args=processargs(defaults,varargin);

      id=obj.checkComposition(ms,'debug',args.debug);
      if ~isempty(args.ref)
        refid=obj.checkComposition(args.ref,'debug',args.debug);
      end
      
      ic=nan(length(id),1);
      relic=ic;
      refic=ic;
      
      p={};
      for i=1:length(id)
        [p{i},r(i),c(i)]=obj.getposition(i);
        if ~isempty(id(i).ic)
          ic(i)=nansum(id(i).ic);
          if ~isempty(args.ref)
            refic(i)=nansum(refid(i).ic);
          else
            relic(i)=id(i).relic;
          end
        end
      end
      if ~isempty(args.ref)
        ti=[ms.name,' vs ',args.ref.name];
      else
        ti=ms.name;
      end
      setfig(ti); clf;

      if ~isempty(args.ref)
        relic=ic./refic;
        relic(refic==0)=nan;

        subplot(333);  
        loglog(refic,ic,'.');
        hold on;
        ax=axis;
        x=logspace(log10(ax(1)),log10(ax(2)));
        plot(x,x+2*sqrt(x),':r');
        plot(x,x-2*sqrt(x),':r');
        xlabel(args.ref.name,'Interpreter','none');
        ylabel(ms.name,'Interpreter','none');
        title('Abs. Ion Count Compare');
        
        subplot(334);
        semilogy([id.mztarget],relic,'.');
        xlabel('m/z');
        ylabel('IC/Ref IC');
        title('Rel IC vs m/z');
        
        subplot(337);
        time=arrayfun(@(z) z.time(1), id);
        semilogy(time,relic,'.');
        xlabel('Elution Time (s)');
        ylabel('IC/Ref IC');
        title('Rel IC vs Time');
        
        subplot(332);
        reftime=arrayfun(@(z) z.time(1), refid);
        h=plot(reftime,time-reftime,'.r');
        hold on;
        sel=relic>args.thresh;
        h(2)=plot(reftime(sel),time(sel)-reftime(sel),'.g');
        ax=axis;
        plot(ax(1:2),-obj.TIMEFUZZ*[1,1],':r');
        plot(ax(1:2),obj.TIMEFUZZ*[1,1],':r');
        xlabel(args.ref.name,'Interpreter','none');
        ylabel('Time diff (s)','Interpreter','none');
        legend(h,{sprintf('RelIC<=%1.g',args.thresh),sprintf('RelIC>%.1g',args.thresh)},'location','best');
        title('Elution Time Compare');
        
      end
      
      fprintf('%s: Located %d compounds with relative ion count >%.2f, %d with >0, out of %d with known elution time, %d total compounds\n', ms.name, sum(ic>=args.thresh), args.thresh, sum(ic>0), sum(isfinite(nanmean(obj.time,2))), length(ic));
      up=unique(p,'sorted');
      
      subplot(331);
      histogram(log10(relic),50)
      ax=axis;
      hold on;
      plot(log10(args.thresh)*[1,1],ax(3:4),'r:');
      if isempty(args.ref)
        title('Ion Count');
        xlabel('log10(Ion Count)');
      else
        title('Relative Ion Count');
        xlabel('log10(Rel. Ion Count)');
      end

      mat=nan(9*4,10*3);
      matic=nan(12,8,10);
      for i=1:length(up)
        row=floor((i-1)/3);
        col=i-row*3-1;
        
        for j=1:max(r)
          for k=min(c):max(c)
            sel=strcmp(p,up{i})&r==j&c==k;
            if any(sel)
              mat(j+row*9+1,(k-1)+col*10)=relic(sel);
              matic(i,j,k)=relic(sel);
            end
          end
        end
      end
      
      subplot(3,3,[5,6,8,9]);
      mat(end+1,:)=nan;
      mat(:,end+1)=nan;
      pcolor(log10(mat));
      shading flat
      axis ij;
      set(gca,'XTick',(1:30)+0.5);
      set(gca,'XTickLabel',arrayfun(@(z) sprintf('%.0f',z), repmat(min(c):max(c),1,3),'UniformOutput',false));
      set(gca,'YTick',(1:36)+0.5);
      lbls=arrayfun(@(z) sprintf('%c',z+'A'-2), repmat(1:max(r)+1,1,4),'UniformOutput',false);
      for i=1:9:length(lbls)
        lbls{i}=' ';
      end
      set(gca,'YTickLabel',lbls);
      set(gca,'TickLength',[0,0])
      hold on;
      for i=1:31
        for j=2:9:37
          plot(i*[1,1],[j,j+8],'-k','LineWidth',1);
        end
      end
      for i=1:10:31
        plot(i*[1,1],[1,37],'-k','LineWidth',4);
      end
      for i=1:37
        plot([1,31],i*[1,1],'-k','LineWidth',1);
      end
      for i=1:9:37
        plot([1,31],i*[1,1],'-k','LineWidth',4);
      end
      for i=1:length(up)
        row=floor((i-1)/3);
        col=i-row*3-1;
        text(col*10+6,row*9+1.55,up{i},'HorizontalAlignment','center','VerticalAlignment','middle');
      end
      
      caxis(log10([args.thresh,nanmax(relic(isfinite(relic(:))))]));
      colorbar;
      if isempty(args.ref)
        title('log10(IC)');
      else
        title('log10(Rel IC)');
      end
      
      if ~isempty(args.ref)
        h=suptitle(sprintf('%s (Ref:%s)',ms.name,args.ref.name));
      else
        h=suptitle(ms.name);
      end
      set(h,'Interpreter','none');
    end
      
    function addFromSDF(obj,ms,sdf,varargin)
    % Add the unique peaks for compounds in given SDF from a particular M/S run
    % Use prior analyses to figure out the expected elution time for each compound
    % or scan all elution times if the no prior data (keep only if a unique peak is determined)
      defaults=struct('debug',false,'group','','contains',[]);  
      args=processargs(defaults,varargin);

      fprintf('Adding data from %s\n',ms.name);
      findex=obj.lookupMS(ms);
      if ~isempty(args.group)
        obj.group{findex}=args.group;
      end
      if isempty(args.contains)
        args.contains=true(size(sdf.sdf));
      end
      % Add all these compounds and mark which ones this file contains
      nindex=[];
      for i=1:length(sdf.sdf)
        s=sdf.sdf(i);
        name=[s.BATCH_PLATE,'-',s.BATCH_WELL];
        mztarget=s.MonoisotopicMass+1;
        nindices(i)=obj.lookupName(name,mztarget,s);
        nindex=nindices(i);
        obj.contains(nindex,findex)=args.contains(i)~=0;   % Mark it as expected to contain
      end
      % Attempt to locate each one uniquely
      for i=1:length(nindices)
        nindex=nindices(i);
        meantime=nanmean(obj.time(nindex,obj.contains(nindex,:)));
        mztarget=obj.mztarget(nindex);
        % Check if the M/Z is unique over the compounds in this file
        samemz=find(obj.contains(:,findex) & abs(obj.mztarget-mztarget)<=obj.MZFUZZ);   % Compounds with same M/Z
        nisomers=1;  nunique=1; ignoreelutetimes=[];
        for i=1:length(samemz)
          if samemz(i)~=nindex
            nisomers=nisomers+1;
            if ~any(isfinite(obj.time(samemz(i),:)))
              % This isomer has unknown elution time
              nunique=nunique+1;
            else
              ignoreelutetimes(end+1)=nanmean(obj.time(samemz(i),:));
            end
          end
        end
        obj.nisomers(nindex,findex)=nisomers;
        obj.nunique(nindex,findex)=nunique;
        
        obj.multihits{nindex,findex}=[];
        obj.numhits(nindex,findex)=0;
        obj.mz(nindex,findex)=nan;
        obj.time(nindex,findex)=nan;
        obj.ic(nindex,findex)=nan;
        if isfinite(meantime)
          id=ms.findcompound(mztarget,'sdf',s,'elutetime',meantime,'timetol',obj.TIMEFUZZ,'mztol',obj.MZFUZZ,'debug',args.debug);
        elseif ~obj.contains(nindex,findex)
          continue;
        elseif nisomers==1
          id=ms.findcompound(mztarget,'sdf',s,'mztol',obj.MZFUZZ,'timetol',obj.TIMEFUZZ,'debug',args.debug);
        elseif nunique==1
          % Have multiple isomers but the other ones all have elution times already
          % Could look for addition peaks
          id=ms.findcompound(mztarget,'sdf',s,'mztol',obj.MZFUZZ,'timetol',obj.TIMEFUZZ,'debug',args.debug,'ignoreelutetimes',ignoreelutetimes);
        else
          if args.debug
            fprintf('Have %d indistinguishable isomers for %s in %s\n', nisomers, obj.names{nindex}, obj.files{findex});
          end
          continue;
        end
        obj.multihits{nindex,findex}=id;
        obj.numhits(nindex,findex)=length(id.mz);
        if length(id.mz)==1
          obj.mz(nindex,findex)=id.mz(1);
          obj.time(nindex,findex)=id.time(1);
          obj.ic(nindex,findex)=id.ic(1);
        elseif length(id.mz)==0
          obj.ic(nindex,findex)=0;
        else
          ;
        end
      end
      minic=1000;
      fprintf('Located %d/%d possible/%d expected compounds and %d unexpected ones with IC>=%.0f\n', sum(obj.contains(:,findex) & isfinite(obj.mz(:,findex))), sum(obj.contains(:,findex) & obj.nunique(:,findex)==1), sum(obj.contains(:,findex)),sum(~obj.contains(:,findex)&isfinite(obj.mz(:,findex))&obj.ic(:,findex)>=minic),minic);
    end
    
    function summary(obj)
    % Summarize data available
      fprintf('Contains %d files, %d compounds (%d with elution time)\n', length(obj.files), length(obj.names), sum(any(isfinite(obj.time'))));
      for i=1:length(obj.files)
        fprintf('%2d %-20.20s %3d/%3d/%3d compounds identified\n', i, obj.files{i}, sum(isfinite(obj.mz(:,i))),sum(obj.nunique(:,i)==1),sum(obj.contains(:,i)))
      end
    end
    
    function x=report(obj)
    % Build table on data by compound
    % Each row is a single compound
    % Data by group
      ugroups=unique(obj.group,'sorted');
      x=[];
      [~,ord]=sort(obj.names);
      for ii=1:length(obj.names)
        i=ord(ii);
        x(ii).name=obj.names{i};
        x(ii).mztarget=obj.mztarget(i);
        
        for j=1:length(ugroups)
          files=find(strcmp(obj.group,ugroups{j})& obj.contains(i,:));
          if ~isempty(files)
          x(ii).mzoffset(j)=nanmean(obj.mz(i,files)')-obj.mztarget(i);
          x(ii).elution(j)=nanmean(obj.time(i,files)');
          x(ii).ioncount(j)=nanmean(obj.ic(i,files)');
          x(ii).nisomers(j)=nanmin(obj.nisomers(i,files)');
          x(ii).nunique(j)=nanmin(obj.nunique(i,files)');
          x(ii).numhits(j)=nanmin(obj.numhits(i,files)');
          f='';
          for k=1:length(files)
            f=[f,obj.files{files(k)},','];
          end
          f=f(1:end-1);  % Remove trailing comma
          x(ii).files{j}=f;
          end
        end
      end
      x=struct2table(x);
    end
    
    function pcolorplot(obj)
    % Heat map of ioncount in matrix of compounds * files
      data=obj.ic;
      % Normalize for xxth percentage
      for i=1:size(data,2)
        data(:,i)=data(:,i)/nanmean(data(:,i));
      end
      data(data<.01)=.01;
      data(data>10)=10;
      setfig('compounds');clf;
      data(end+1,:)=nan;
      data(:,end+1)=nan;
      pcolor(log10(data)');
      shading flat;
      colorbar;
      xlabel('Compound');
      ylabel('File');
      set(gca,'YTick',(1:length(obj.files))+0.5);
      set(gca,'YTickLabels',strrep(obj.files,'.mzXML',''));
      set(gca,'ticklabelinterpreter','none');
    end
    
    function scaling=getscaling(obj,f1,f2)
    % Get ioncounts in f2 relative to f1
      ratio=obj.ic(:,f2)./obj.ic(:,f1);
      scaling=nanmedian(ratio(ratio>0));
    end
    
    function checkmzoffset(obj)
    % Check whether mzoffset used when reading mass spec files should be changed
      fprintf('File          Additional Offset\n');
      ugroups=unique(obj.group);
      for j=1:length(ugroups)
        all=[];
        for i=1:length(obj.files)
          if strcmp(obj.group{i},ugroups{j})
            err=nanmedian(obj.mz(:,i)-obj.mztarget);
            fprintf('%-20.20s  %8.4f\n', obj.files{i}, err);
            all(end+1)=err;
          end
        end
        sel=strcmp(obj.group,ugroups{j});
        fprintf('%-20.20s  %8.4f over group\n', ugroups{j}, nanmedian(all));
      end
    end
    
    function checktimeoffset(obj,sel)
    % Check time alignments
      meantime=nanmean(obj.time,2);
      [meantime,ord]=sort(meantime);
      time=obj.time(ord,:);
      reltime=nan(size(time));
      for i=1:size(time,2)
        tother=time;  tother(:,i)=nan;
        reltime(:,i)=time(:,i)-nanmean(tother,2);
      end
      setfig('timeoffset');clf;
      tiledlayout('flow');
      for i=1:length(sel)
        nexttile
        plot(meantime,reltime(:,sel(i)),'.');
        use=isfinite(meantime) & isfinite(reltime(:,sel(i)));
        p=polyfit(meantime(use),reltime(use,sel(i)),1)
        hold on;
        plot(meantime,polyval(p,meantime),'-r');
        xlabel('Mean time');
        ylabel('Individual run times - meantime');
        title(sprintf('%s m=%.4f, b=%.4f',obj.files{sel(i)},p));
      end
    end
    
    
    function plotcompare(obj,f1,f2)
    % Plot comparison of each compound the occurs in both f1 and f2
      ti=sprintf('%s vs %s',obj.files{f1}, obj.files{f2});
      setfig(ti);clf;
      subplot(221)
      plot(obj.mz(:,f1)-obj.mztarget,obj.mz(:,f2)-obj.mztarget,'o');
      hold on;
      ax=axis;
      plot(ax(1:2),obj.MZFUZZ*[1,1],'r:');
      plot(ax(1:2),-obj.MZFUZZ*[1,1],'r:');
      ax=axis;
      plot(obj.MZFUZZ*[1,1],ax(3:4),'r:');
      plot(-obj.MZFUZZ*[1,1],ax(3:4),'r:');
      xlabel(obj.files{f1},'Interpreter','none');
      ylabel(obj.files{f2},'Interpreter','none');
      title('m/z offset');

      subplot(222);
      plot(mean(obj.time(:,[f1,f2]),2),diff(obj.time(:,[f1,f2]),[],2),'o');
      hold on;
      ax=axis;
      plot(ax(1:2),obj.TIMEFUZZ*[1,1],'r:');
      plot(ax(1:2),-obj.TIMEFUZZ*[1,1],'r:');
      xlabel('Mean (s)');
      ylabel('Diff (s)');
      title('Elution Times');
      
      ratio=obj.ic(:,f2)./obj.ic(:,f1);
      scaling=nanmedian(ratio(ratio>0));
      ratio=ratio/scaling;
      ratio(ratio==0)=.001;
      ratio(obj.ic(:,f1)==0 & obj.ic(:,f2)~=0)=1000;

      subplot(223);
      x=obj.ic(:,f1); y=obj.ic(:,f2);
      x(x==0)=1;
      y(y==0)=1;
      loglog(x,y,'o');
      hold on;
      ax=axis;
      plot(ax(1:2),ax(1:2)*scaling,':');
      xlabel(obj.files{f1},'Interpreter','none');
      ylabel(obj.files{f2},'Interpreter','none');
      title('Ion Counts');
      
      subplot(224);
      %      semilogy(obj.mztarget,max(ratio,.01),'o');
      histogram(log10(ratio),20);
      xlabel('log10(IC2/IC1)');
      ylabel('N');
      title(sprintf('Ion Count Ratio (scaling=%.2f)',scaling));
      
      h=suptitle(sprintf('%s N=%d',ti,sum(all(isfinite(obj.ic(:,[f1,f2])),2))));
      set(h,'Interpreter','none');
    end

    function plotmap(obj)
    % Plot map of compounds in m/z vs time space
      setfig('Compound Map');clf;
      etime=[];
      for i=1:length(obj.mztarget)
        etime(i,1)=nanmin(obj.time(i,:));
        etime(i,2)=nanmax(obj.time(i,:));
        etime(i,3)=nanmean(obj.time(i,:));
        mz(i,1)=nanmin(obj.mz(i,:));
        mz(i,2)=nanmax(obj.mz(i,:));
      end
      plot(obj.mztarget,etime(:,3),'.b');
      hold on;
      ngood=0;
      for i=1:size(etime,1)
        if isnan(etime(i,3))
          continue;
        end
        overlap=sum(etime(:,1)<etime(i,2) & etime(:,2)>etime(:,1) & mz(:,1)<mz(i,2) & mz(:,2)>mz(i,1));
        if overlap>1
          col='r';
        else
          col='g';
          ngood=ngood+1;
        end
        plot(mz(i,[1,2,2,1,1]),etime(i,[1,1,2,2,1]),col);
      end
      xlabel('m/z');
      ylabel('Elution time (s)');
      title(sprintf('Compound Map (%d distinguishable/%d identified)',ngood,sum(isfinite(etime(:,3)))));
    end
    
    function getinfo(obj,name,varargin)
      defaults=struct('mzdata',[]);
      args=processargs(defaults,varargin);

      if ischar(name)
        ind = obj.find(name);
      else
        ind=name;
      end
      meanic=nanmean(obj.ic(ind,obj.contains(ind,:)));
      meant=nanmean(obj.time(ind,obj.contains(ind,:)));
      fprintf('%s (%d): m/z=%8.4f t=%7.2f meanic=%.0f\n',obj.shortnames{ind},ind, obj.mztarget(ind),meant,meanic);
      if ~isempty(args.mzdata)
        setfig(obj.shortnames{ind});
        t=tiledlayout('flow');
        title(t,sprintf('%s m/z=%.4f t=%.0f',obj.shortnames{ind},obj.mztarget(ind),meant));
      end
      
      for j=1:length(obj.files)
        if ~obj.contains(ind,j)
          % TODO: Could list false positives here
          continue;
        end
        fprintf('%-15.15s: nisomers=%d, nunique=%d, nhits=%-2d ',obj.files{j},obj.nisomers(ind,j), obj.nunique(ind,j), obj.numhits(ind,j));
        m=obj.multihits{ind,j};
        if ~isempty(m)
          for k=1:length(m.mz)
            fprintf('[mz=%8.4f, t=%.2f, ic=%.0f] ', m.mz(k), m.time(k), m.ic(k));
          end
        end
        if ~isempty(args.mzdata)
          [ic,mz,t]=args.mzdata{j}.mzscan(obj.mztarget(ind),'mztol',obj.MZFUZZ);
          nexttile;
          plot(t,ic);
          ax=axis;
          hold on;
          plot(meant+obj.TIMEFUZZ*[1,1],ax(3:4),':b');
          plot(meant-obj.TIMEFUZZ*[1,1],ax(3:4),':b');
          ylabel('Ion Count');
          yyaxis right
          plot(t,mz,'r');
          hold on;
          ax=axis;
          axis(ax);
          plot(ax(1:2),obj.mztarget(ind)+obj.MZFUZZ*[1,1],':r');
          plot(ax(1:2),obj.mztarget(ind)-obj.MZFUZZ*[1,1],':r');
          ylabel('M/Z');
          if isfinite(meant)
            ax(1:2)=meant+obj.TIMEFUZZ*2*[-1,1];
          end
          ax(3:4)=obj.mztarget(ind)+obj.MZFUZZ*2*[-1,1];
          axis(ax);
          title(sprintf('%s',args.mzdata{j}.name));
        end
        fprintf('\n');
      end
      if isfinite(meanic) && isfinite(meant)
        % False positives
        minic=meanic/10;
        fprintf('False positives with  m/z in [%.3f,%.3f], T in [%.0f,%.0f], IC >= %.0f:\n',...
                obj.mztarget(ind)+obj.MZFUZZ*[-1,1],...
                nanmean(obj.time(ind,obj.contains(ind,:)))+obj.TIMEFUZZ*[-1,1],...
                minic);
        for j=1:length(obj.files)
          if ~obj.contains(ind,j) && obj.ic(ind,j)>=minic
            fprintf('%-15.15s: mz=%8.4f, t=%7.2f, ic=%5.0f\n',obj.files{j},obj.mz(ind,j), obj.time(ind,j), obj.ic(ind,j));
          end
        end
      end
    end
    
  end
end
