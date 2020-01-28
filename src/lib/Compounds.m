% Data structure to hold information about compounds located in mass spec runs
classdef Compounds < handle
  properties
    names;   % names{i} - Name of compound i (e.g.'51B07')
    adduct;  % adduct(i) - Name of adduct of compound (e.g. 'H')
    sdf;     % sdf{i} - SDF data for compound i
    mztarget; % mztarget(i) - target m/z for compound i with adduct(i)
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
    ADDUCTS=struct('name',{'H','Na'},'mass',{1.007825035,22.9897677});
  end
  
  methods
    function obj=Compounds()
      obj.multihits={};
      obj.contains=false(0,0);
    end

    function t=meantime(obj,i)
      if nargin>1
        t=nanmean(obj.time(i,obj.contains(i,:)));
      else
        t=[];
        for i=1:size(obj.time,1)
          t(i)=nanmean(obj.time(i,obj.contains(i,:)));
        end
        t=t';
      end
    end
    
    function ind=find(obj,name,adduct)
    % Find given name
    % Try both format stored in names (e.g. 'CDIV0051-B07') and short format (e.g. '51B07')
      if nargin<3
        c=strsplit(name,'+');
        if length(c)==2
          name=c{1};
          adduct=c{2};
        else
          adduct='H';
        end
      end
      ind=find(strcmp(name,obj.names) & stcmp(obj.adduct,adduct));
      if length(ind)~=1
        error('Compound %s+%s not found',name,adduct);
      end
    end
      
    function nindex=lookupName(obj,name, mztarget)
    % Find the index of compound by name, create if missing
      if ismember(name,obj.names)
        nindex=find(strcmp(name,obj.names));
      else
        obj.names{end+1,1}=name;
        nindex=length(obj.names);
        if nargin>=3
          obj.mztarget(nindex,1)=mztarget;
        else
          obj.mztarget(nindex,1)=nan;
        end
        obj.mz(nindex,:)=nan;
        obj.time(nindex,:)=nan;
        obj.ic(nindex,:)=nan;
        obj.contains(nindex,:)=false;
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
          id(i)=ms.findcompound(mztarget,'elutetime',meantime(i),'timetol',obj.TIMEFUZZ,'mztol',obj.MZFUZZ,'debug',args.debug);
          if isempty(id(i).ic)
            id(i).ic=0;
            id(i).time=nan;
            id(i).mz=nan;
          end
          if length(isomers)>1 && sum(id(i).ic)>0
            fprintf('Isomer at compounds %s with ion count = %.0f\n', sprintf('%d,',isomers),sum(id(i).ic));
            id(i).ic=nan;   % Can't use it
          end
        end
      end
      
      for i=1:length(obj.mztarget)
        id(i).name=obj.names{i};
        id(i).adduct=obj.adduct{i};
        id(i).relic=sum(id(i).ic)/nanmax(obj.ic(i,:));
      end
    end
    
    function [matic,id,refid]=plotComposition(obj,ms,varargin)
      defaults=struct('debug',false,'thresh',0.02,'ref',[]);  
      args=processargs(defaults,varargin);

      id=obj.checkComposition(ms,'debug',args.debug);
      if ~isempty(args.ref)
        refid=obj.checkComposition(args.ref,'debug',args.debug);
      else
        refid=nan;
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
              matic(i,j,k-1)=relic(sel);
            end
          end
        end
      end
      
      subplot(3,3,[5,6,8,9]);
      data=mat;
      data(end+1,:)=nan;
      data(:,end+1)=nan;
      pcolor(log10(data));
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
      
    function addCompound(obj,name,mass,adduct)
    % Add a compound
      if any(strcmp(obj.names,name) & strcmp(obj.adduct,adduct))
        error('%s+%s already added', name, adduct);
      end
      addsel=find(strcmp({obj.ADDUCTS.name},adduct));
      if isempty(addsel)
        error('Unknown adduct: %s',adduct);
      end

      obj.names{end+1}=name;
      nindex=length(obj.names);
      obj.adduct{nindex}=adduct;
      obj.mztarget{nindex}=mass+obj.ADDUCTS(addsel).mass;
      obj.mz(nindex,:)=nan;
      obj.time(nindex,:)=nan;
      obj.ic(nindex,:)=nan;
      obj.contains(nindex,:)=false;
    end
    
    function addCompoundsFromSDF(obj,sdf,adduct)
    % Add all the compounds in the given SDF file using a name formed from the PLATE and WELL
      if nargin<3
        adduct='H';
      end
      for i=1:length(sdf.sdf)
        s=sdf.sdf(i);
        name=sprintf('%d%s',str2num(s.BATCH_PLATE(5:end)),s.BATCH_WELL);
        obj.addCompound(name,s.MonoisotopicMass,adduct);
      end
    end
    
    function addMS(obj,ms,varargin)
    % Add the unique peaks for compounds in given SDF from a particular M/S run
    % Use prior analyses to figure out the expected elution time for each compound
    % or scan all elution times if the no prior data (keep only if a unique peak is determined)
      defaults=struct('debug',false,'group','','contains',{{}},'mztol',[],'timetol',[]);  
      args=processargs(defaults,varargin);

      if isempty(args.mztol)
        args.mztol=obj.MZFUZZ;
      end
      if isempty(args.timetol)
        args.timetol=obj.TIMEFUZZ;
      end
      
      fprintf('Adding data from %s\n',ms.name);
      findex=obj.lookupMS(ms);
      if ~isempty(args.group)
        obj.group{findex}=args.group;
      end
      if isempty(args.contains)
        obj.contains(:,findex)=true;
      else
        badnames=setdiff(args.contains,obj.names);
        if ~isempty(badnames)
          error('contains list has nonexistent compounds names: ',strjoin(badnames,','));
        end
        obj.contains(:,findex)=ismember(obj.names,args.contains);
      end

      % Attempt to locate each one uniquely
      for nindex=1:length(obj.mztarget)
        meantime=nanmean(obj.time(nindex,obj.contains(nindex,:)));
        mztarget=obj.mztarget(nindex);
        % Check if the M/Z is unique over the compounds in this file
        samemz=find(obj.contains(:,findex) & abs(obj.mztarget-mztarget)<=args.mztol);   % Compounds with same M/Z
        nisomers=1;  nunique=1; ignoreelutetimes=[]; nident=1;
        for i=1:length(samemz)
          if samemz(i)~=nindex
            nisomers=nisomers+1;
            if ~any(isfinite(obj.time(samemz(i),:)))
              % This isomer has unknown elution time
              nunique=nunique+1;
            else
              if abs(nanmean(obj.time(samemz(i),:)) - meantime) < args.timetol
                if obj.contains(nindex,findex) && nindex<samemz(i)
                  fprintf('%s and %s have same m/z, same elution time\n', obj.names{nindex}, obj.names{samemz(i)});
                end
                nunique=nunique+1;  % Same elution times
                nident=nident+1;
              else
                ignoreelutetimes(end+1)=nanmean(obj.time(samemz(i),:));
              end
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
        if isfinite(meantime) && nident==1
          id=ms.findcompound(mztarget,'elutetime',meantime,'timetol',args.timetol,'mztol',args.mztol,'debug',args.debug);
        elseif ~obj.contains(nindex,findex)
          continue;
        elseif nisomers==1
          id=ms.findcompound(mztarget,'mztol',args.mztol,'timetol',args.timetol,'debug',args.debug);
        elseif nunique==1
          % Have multiple isomers but the other ones all have elution times already
          % Could look for addition peaks
          id=ms.findcompound(mztarget,'mztol',args.mztol,'timetol',args.timetol,'debug',args.debug,'ignoreelutetimes',ignoreelutetimes);
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
        fprintf('%2d %-20.20s %3d/%3d/%3d compounds identified/unique/total\n', i, obj.files{i}, sum(obj.contains(:,i) & isfinite(obj.mz(:,i))),sum(obj.contains(:,i) & obj.nunique(:,i)==1),sum(obj.contains(:,i)))
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
        x(ii).adduct=obj.adduct{i};
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
            fit=robustfit(obj.mztarget,obj.mz(:,i));
            fprintf('%-20.20s  %8.4f [%8.4f@%3.0f - %8.4f@%3.0f]\n', obj.files{i}, err,fit(1)+(fit(2)-1)*min(obj.mztarget),min(obj.mztarget),fit(1)+(fit(2)-1)*max(obj.mztarget),max(obj.mztarget));
            all(end+1)=err;
          end
        end
        sel=strcmp(obj.group,ugroups{j});
        fit=robustfit(obj.mztarget,nanmean(obj.mz(:,sel),2));
        fprintf('%-20.20s  %8.4f [%8.4f@%3.0f - %8.4f@%3.0f] over group\n', ugroups{j}, nanmedian(all), fit(1)+(fit(2)-1)*min(obj.mztarget),min(obj.mztarget),fit(1)+(fit(2)-1)*max(obj.mztarget),max(obj.mztarget));
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
    
    function map=checktimeoffset2(obj,obj2,sel1,sel2)
    % Compare time offsets between two structures with same compound list
      if nargin<3
        sel1=1;
      end
      if nargin<4
        sel2=1;
      end
      assert(all(obj.mztarget==obj2.mztarget));
      setfig('checktimeoffset2');clf;
      subplot(211);
      t1=obj.time(:,sel1); t2=obj2.time(:,sel2);
      plot(t1,t2,'o'); hold on;
      map=piecewise(t1,t2,obj.TIMEFUZZ,4);
      pred=interp1(map(:,1),map(:,2),t1);
      plot(t1,pred,'r-');
      subplot(212);
      resid=t2-pred;
      plot(t1,resid,'o');
      ylabel('Residual');
      rmse=16;
      for i=1:4
        outlier=abs(resid)>rmse*4;
        rmse=sqrt(nanmean(resid(~outlier).^2));
      end
      hold on; plot(t1(outlier),resid(outlier),'or');
      fprintf('RMSE(resid) = %.1f (outliers are >=%.0f)\n',rmse, min([inf;abs(resid(outlier))])); 
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
      fprintf('%s (%d): m/z=%8.4f t=%7.2f meanic=%.0f\n',obj.names{ind},ind, obj.mztarget(ind),meant,meanic);
      if ~isempty(args.mzdata)
        setfig(obj.names{ind});
        t=tiledlayout('flow');
        title(t,sprintf('%s m/z=%.4f t=%.0f',obj.names{ind},obj.mztarget(ind),meant));
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
          nexttile;
          obj.plotscan(ind,args.mzdata{j});
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
            if ~isempty(args.mzdata)
              nexttile;
              obj.plotscan(ind,args.mzdata{j});
            end
          end
        end
      end
    end

    function plotscan(obj,ind,mzdata)
      meanic=nanmean(obj.ic(ind,obj.contains(ind,:)));
      meant=nanmean(obj.time(ind,obj.contains(ind,:)));
      [ic,mz,t]=mzdata.mzscan(obj.mztarget(ind),'mztol',obj.MZFUZZ);
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
      title(sprintf('%s',mzdata.name));
    end
    
    function listunexpected(obj,varargin)
    % List peaks that shouldn't be present, in order of descending ion count
      defaults=struct('nlist',20);
      args=processargs(defaults,varargin);

      ic=obj.ic;
      ic(obj.contains)=0;
      [~,ord]=sort(ic(:),'desc','MissingPlacement','last');
      for i=1:args.nlist
        [ii,ij]=ind2sub(size(obj.contains),ord(i));
        fprintf('%12.12s:%-7.7s IC=%8.0f', obj.files{ij}, obj.names{ii}, ic(ord(i)));
        alias=find(obj.contains(:,ij) & abs(obj.mztarget-obj.mztarget(ii))<=obj.MZFUZZ & abs(obj.meantime-obj.meantime(ii))<=obj.TIMEFUZZ);
        if ~isempty(alias)
          fprintf(' Indistinguishable from ');
          for j=1:length(alias)
            fprintf(' %s',obj.names{alias(j)});
          end
        end
        fprintf('\n');
      end
    end
  end
end
