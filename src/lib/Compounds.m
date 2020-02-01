% Data structure to hold information about compounds located in mass spec runs
classdef Compounds < handle
  properties
    names;   % names{i} - Name of compound i (e.g.'51B07')
    adduct;  % adduct(i) - Name of adduct of compound (e.g. 'H')
    mztarget; % mztarget(i) - target m/z for compound i with adduct(i)
    files;   % files{j} - Mass spec filename j
    moles;    % moles(j) - Moles of each compound loaded in run j
    group;   % group{j} - name of group this file belongs to
    contains;  % contains(i,j) is true if we expected to find compound i in file j
             % For compounds uniquely located in a particular run:
    mz;    % mz(i,j) contains the observed m/z peak on compound i, from file j
    time;    % time(i,j) contains the elution time on compound i, from file j
    filetime;    % filetime(i,j) contains the elution time on compound i, from file j without any time remapping
    ic;    % ic(i,j) contains the total ion count for compound i, from file j
    normic;	% Normalized ion count (by file and by target)
    multihits;
    tsens;    % tsens(i) is the relative sensitivity to target i
    fsens;    % fsens(j) is the relative sensitivity for file j
  end
  
  properties(Constant)
    MZFUZZ=0.004; % 0.003;
    TIMEFUZZ=50;   % in seconds
    ADDUCTS=struct('name',{'M+H','M+Na'},'mass',{1.007825035,22.9897677});
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
          adduct='';
        end
      end
      if isempty(adduct) 
        ind=find(strcmp(name,obj.names));
        %if length(ind)>1
        % fprintf('Warning: Compound %s has multiple adducts: %s',name,strjoin(obj.adduct(ind),','));
        %end
      else
        ind=find(strcmp(name,obj.names) & strcmp(obj.adduct,adduct));
      end
      if length(ind)==0
        error('Compound %s+%s not found',name,adduct);
      end
    end
      
    function findex=lookupMS(obj,ms)
    % Find the index of a file by name, create if missing
      if ismember(ms.path,obj.files)
        findex=find(strcmp(ms.path,obj.files));
      else
        obj.files{1,end+1}=ms.path;
        findex=length(obj.files);
        obj.mz(:,findex)=nan;
        obj.time(:,findex)=nan;
        obj.filetime(:,findex)=nan;
        obj.ic(:,findex)=nan;
        obj.contains(:,findex)=false;
        obj.moles(findex)=ms.moles;
      end
    end

    function [p,r,c]=getposition(obj,i)
    % Get plate (str),row (int),col(int) of compound i
      A=sscanf(obj.names{i},'%d%c%2d');
      p=A(1);
      r=A(2)-'A'+1;
      c=A(3);
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
      
      p=[];r=[];c=[];
      for i=1:length(id)
        [p(i),r(i),c(i)]=obj.getposition(i);
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
            sel=p==up(i)&r==j&c==k;
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
        text(col*10+6,row*9+1.55,sprintf('%d',up(i)),'HorizontalAlignment','center','VerticalAlignment','middle');
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
      obj.mztarget(nindex)=mass+obj.ADDUCTS(addsel).mass;
      if size(obj.mz,2)>0
        obj.mz(nindex,:)=nan;
        obj.time(nindex,:)=nan;
        obj.filetime(nindex,:)=nan;
        obj.ic(nindex,:)=nan;
        obj.contains(nindex,:)=false;
        obj.multihits{nindex,:}=[];
      else
        obj.mz=nan(nindex,0);
        obj.time=nan(nindex,0);
        obj.filetime=nan(nindex,0);
        obj.ic=nan(nindex,0);
        obj.contains=false(nindex,0);
        obj.multihits=cell(nindex,0);
      end
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
    
    function assignTimes(obj,varargin)
    % Collect all the MS data and assign elution times to compounds
    % For each compounds, 
    %   build list of observations across massspec files (elution time,whether compound is expected)
    %   find elution time with highest correlation of hit vs. expected
      defaults=struct('debug',false,'timetol',obj.TIMEFUZZ,'minhits',4);
      args=processargs(defaults,varargin);

      fprintf('assignTimes:\n');
      obj.ic=nan(length(obj.names),length(obj.files));
      obj.mz=nan(length(obj.names),length(obj.files));
      obj.time=nan(length(obj.names),length(obj.files));
      obj.filetime=nan(length(obj.names),length(obj.files));

      for i=1:length(obj.names)
        fprintf('%s: ',obj.names{i});
        etimes=[];ic=[];cont=false(0,0);
        for j=1:length(obj.files)
          m=obj.multihits{i,j};
          if ~isempty(m)
            etimes=[etimes,m.time];
            ic=[ic,m.ic];
            cont=[cont,repmat(obj.contains(i,j),size(m.ic))];
          end
        end
        contains=obj.contains(i,:);
        if length(etimes)>1
          t=clusterdata(etimes','Criterion','distance','cutoff',args.timetol,'linkage','centroid');
          meantime=[];
          for j=1:max(t)
            meantime(j)=mean(etimes(t==j));
          end
        
          hits=false(max(t),length(obj.files));
          for j=1:length(obj.files)
            m=obj.multihits{i,j};
            if ~isempty(m)
              [~,~,ib]=intersect(m.time,etimes);
              hits(t(ib),j)=true;
            end
          end
          matches=[];
          for j=1:size(hits,1)
            matches(j,1)=sum(hits(j,:)&contains);
            matches(j,2)=sum(~hits(j,:)&~contains);
          end
          score=matches*[5;1];
          [~,best]=max(score);
        elseif length(etimes)==1
          best=1;
          meantime=etimes;
        else
          fprintf('no hits\n');
          continue;
        end
        fprintf('%d/%d hits with %d false hits at T=%.0f ', matches(best,1),sum(contains), sum(~contains)-matches(best,2), meantime(best));
        if matches(best,1)<args.minhits
          fprintf('too few hits');
        elseif sum(score==max(score))>1
          fprintf('ambiguous');
        else
          % Construct consensus view
          for j=1:length(obj.files)
            m=obj.multihits{i,j};
            if ~isempty(m)
              sel=abs(m.time-meantime(best))<args.timetol;
              if sum(sel)>0
                obj.ic(i,j)=sum(m.ic(sel));
                obj.mz(i,j)=sum(m.mz(sel).*m.ic(sel))/obj.ic(i,j);
                %obj.filemz(i,j)=sum(m.filemz(sel).*m.ic(sel))/obj.ic(i,j);
                obj.time(i,j)=sum(m.time(sel).*m.ic(sel))/obj.ic(i,j);
                obj.filetime(i,j)=sum(m.filetime(sel).*m.ic(sel))/obj.ic(i,j);
              end
            end
          end
        end
        fprintf('\n');
        
        if args.debug
          setfig('assignTimes');clf;
          semilogy(etimes(~cont),ic(~cont),'or');
          hold on;
          semilogy(etimes(cont),ic(cont),'xb');
          ax=axis;
          for j=1:length(meantime)
            plot((meantime(j)-args.timetol)*[1,1],ax(3:4),':g');
            plot((meantime(j)+args.timetol)*[1,1],ax(3:4),':g');
          end
          plot(meantime(best)*[1,1],ax(3:4),'g');
          xlabel('Elution time');
          ylabel('Ion count');
          title(obj.names{i});
          figure(gcf);
          keyboard;
        end
      end
    end
    
    function x=icPerMole(obj)
    % Get median value of ion count per mole of target loaded
      x=nanmedian(obj.ic,1)./obj.moles;
    end
    
    function addMS(obj,ms,varargin)
    % Add the unique peaks for compounds in given SDF from a particular M/S run
    % Use prior analyses to figure out the expected elution time for each compound
    % or scan all elution times if the no prior data (keep only if a unique peak is determined)
      defaults=struct('debug',false,'group','','contains',{{}},'mztol',[],'timetol',[],'map',[]);
      % mzmap(i,2) - piecewise linear map for M/Z; mzmap(:,1) is true M/Z, mzmap(:,2) is for values in ms file
      % timemap(i,2) - piecewise linear map for elution times; timemap(:,1) is "standard" elution time
      args=processargs(defaults,varargin);

      if isempty(args.mztol)
        args.mztol=obj.MZFUZZ;
      end
      if isempty(args.timetol)
        args.timetol=obj.TIMEFUZZ;
      end
      if isempty(args.map)
        args.map=struct('mz',[0 0; 1 1 ],'time',[0 0; 1 1 ]);
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
      maxic=[];   % Maximum IC per target
      for i=1:length(obj.mztarget)
        mztargetMS=interp1(args.map.mz(:,1),args.map.mz(:,2),obj.mztarget(i),'linear','extrap');
        id=ms.findcompound(mztargetMS,'mztol',args.mztol,'timetol',args.timetol,'debug',args.debug);
        % Map M/Z, time back to global values
        id.filemz=id.mz;
        id.filetime=id.time;
        id.mz=interp1(args.map.mz(:,2),args.map.mz(:,1),id.mz,'linear','extrap');
        id.time=interp1(args.map.time(:,2),args.map.time(:,1),id.time,'linear','extrap');
        obj.multihits{i,findex}=id;
        maxic(i)=max([0,id.ic]);
      end
      ice=maxic(obj.contains(:,findex));
      icu=maxic(~obj.contains(:,findex));
      if isempty(icu) || nanmedian(icu)<200
        minic=prctile(ice(ice>0),10);
      else
        minic=sqrt(median(ice(ice>0))*median(icu(icu>0)));
      end
      p=[25,50,75];
      fprintf('Expected have IC=[%s], unexpected have IC=[%s] (%s) @%.0f: %.1f%%,%.1f%% \n', sprintf('%.0f ',prctile(ice,p)), sprintf('%.0f ',prctile(icu,p)), sprintf('%d%% ',p),minic,100*mean(ice>minic),100*mean(icu>minic));
      nhits=nan(length(obj.mztarget),1);
      for i=1:length(obj.mztarget)
        nhits(i)=sum(obj.multihits{i,findex}.ic>minic);
      end
      fprintf('Have hits for %d/%d (with %d unique) expected compounds and %d unexpected ones with IC>=%.0f\n', sum(obj.contains(:,findex) & nhits>0), sum(obj.contains(:,findex)), sum(obj.contains(:,findex) & nhits==1), sum(~obj.contains(:,findex)&nhits>0),minic);
    end
    
    function summary(obj)
    % Summarize data available
      fprintf('summary():\n');
      fprintf('Contains %d files, %d compounds (%d with elution time)\n', length(obj.files), length(obj.names), sum(any(isfinite(obj.time'))));
      for i=1:length(obj.files)
        [~,filename,~]=fileparts(obj.files{i});
        fprintf('%2d %-20.20s %3d/%3d compounds identified/total\n', i,filename, sum(obj.contains(:,i) & isfinite(obj.mz(:,i))),sum(obj.contains(:,i)));
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
          f='';
          for k=1:length(files)
            [~,filename,~]=fileparts(obj.files{k});
            f=[f,filename,','];
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
        data(:,i)=data(:,i)/prctile(data(:,i),80);
      end
      data(data<.01)=0;
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
      filenames={};
      for i=1:length(obj.files)
        [~,filenames{i},~]=fileparts(obj.files{i});
      end
      set(gca,'YTickLabels',filenames);
      set(gca,'ticklabelinterpreter','none');
    end
    
    function scaling=getscaling(obj,f1,f2)
    % Get ioncounts in f2 relative to f1
      ratio=obj.ic(:,f2)./obj.ic(:,f1);
      scaling=nanmedian(ratio(ratio>0));
    end
    
    function checkmzoffset(obj)
    % Check whether mzoffset used when reading mass spec files should be changed
      fprintf('checkmzoffset()\n');
      dirs={};
      for i=1:length(obj.files)
        dirs{i}=fileparts(obj.files{i});
      end
      udirs=unique(dirs);
      setfig('checkmzoffset');clf;
      tiledlayout('flow');
      rng=[min(obj.mztarget),max(obj.mztarget)];
      fprintf('File     Added Offset(*1e4) [%5.0f, %5.0f]\n',rng);
      for j=1:length(udirs)
        nexttile;
        all=[];
        h=[];
        leg={};
        for i=1:length(obj.files)
          if strcmp(dirs{i},udirs{j})
            err=nanmedian(obj.mz(:,i)'-obj.mztarget);
            fit=robustfit(obj.mztarget,obj.mz(:,i));
            [~,filename]=fileparts(obj.files{i});
            fprintf('%-20.20s  %5.1f [%5.1f, %5.1f]\n',filename, err*1e4,1e4*(fit(1)+(fit(2)-1)*rng));
            h(end+1)=plot(obj.mztarget,10000*(obj.mz(:,i)'-obj.mztarget),'o');
            hold on;
            plot(rng,10000*(fit(1)+(fit(2)-1)*rng),'-','Color',get(h(end),'Color'));
            leg{end+1}=filename;
            all(end+1)=err;
          end
        end
        sel=strcmp(dirs,udirs{j});
        fit=robustfit(obj.mztarget,nanmean(obj.mz(:,sel),2));
        [~,dirname]=fileparts(udirs{j});
        fprintf('%-20.20s  %5.1f [%5.1f, %5.1f] over %s\n', '', 1e4*nanmedian(all), 1e4*(fit(1)+(fit(2)-1)*rng),dirname);
        h(end+1)=plot(rng,10000*(fit(1)+(fit(2)-1)*rng),'k','linewidth',2);
        plot(rng,10000*(fit(1)+(fit(2)-1)*rng+obj.MZFUZZ),'k:');
        plot(rng,10000*(fit(1)+(fit(2)-1)*rng-obj.MZFUZZ),'k:');
        leg{end+1}='All';
        legend(h,leg,'location','best');
        xlabel('True M/Z');
        ylabel('File-True M/Z * 10000');
        title(udirs{j});
      end
    end
    
    function map=checktime(obj,ref,varargin)
    % Check all multiple hits of ref against each other
      defaults=struct('debug',false,'timetol',obj.TIMEFUZZ);
      args=processargs(defaults,varargin);

      fprintf('checktime(%d,''timetol'',%.0f)\n',ref,args.timetol);
      t = nan(size(obj.mz));
      ic = nan(size(obj.mz));
      for i=1:length(obj.mztarget)
        for j=1:length(obj.files)
          if obj.contains(i,j)
            m=obj.multihits{i,j};
            if ~isempty(m) && ~isempty(m.ic)
              [ic(i,j),ord]=max(m.ic);
              t(i,j)=m.filetime(ord);
            end
          end
        end
      end
      dirs={};
      for i=1:length(obj.files)
        dirs{i}=fileparts(obj.files{i});
      end
      udirs=unique(dirs);

      tref=t(:,ref);
      trefs=sort(tref);
      setfig('checktime');clf;
      tiledlayout('flow');
      map=[];
      for j=1:length(udirs);
        nexttile;
        h=[];leg={};
        for i=1:length(obj.files)
          if strcmp(dirs{i},udirs{j})
            t2=t(:,i);
            fit=piecewise(tref,t2,args.timetol,2);
            pred=interp1(fit(:,1),fit(:,2),trefs,'linear','extrap');
            [~,filename]=fileparts(obj.files{i});
            fprintf('%-20.20s  [%s]\n',filename, sprintf('(%4.0f@%4.0f) ',fit'));
            h(end+1)=plot(tref,t2-tref,'o');
            hold on;
            plot(trefs,pred-trefs,'-','Color',get(h(end),'Color'));
            leg{end+1}=filename;
          end
        end
        sel=strcmp(dirs,udirs{j});
        fit=piecewise(tref,nanmedian(t(:,sel),2),args.timetol,2);
        [~,dirname]=fileparts(udirs{j});
        fprintf('%-20.20s  [%s] over %s\n\n', '', sprintf('(%4.0f@%4.0f) ',fit'),dirname);
        pred=interp1(fit(:,1),fit(:,2),trefs,'linear','extrap');
        h(end+1)=plot(trefs,pred-trefs,'k','linewidth',2);
        leg{end+1}='All';
        plot(trefs,pred-trefs+obj.TIMEFUZZ,'k:','linewidth',1);
        plot(trefs,pred-trefs-obj.TIMEFUZZ,'k:','linewidth',1);
        ax=axis;
        ax(3)=min(pred-trefs-obj.TIMEFUZZ*2);
        ax(4)=max(pred-trefs+obj.TIMEFUZZ*2);
        axis(ax);
        legend(h,leg,'location','best');
        xlabel('Reference time');
        ylabel('File-Ref time');
        title(udirs{j});
        map=[map,struct('dir',udirs{j},'map',fit)];
      end
    end
    
    function checksensitivity(obj,ref)
    % Check sensitivity by file relative to ref file
      sens=[]; lbl={};
      for i=1:length(obj.files)
        sel=all(obj.contains(:,[i,ref]),2);
        sens(i)=nanmedian(obj.ic(sel,i)./obj.ic(sel,ref));
        [~,lbl{i}]=fileparts(obj.files{i});
      end
      setfig('File Sensitivity');clf;
      bar(sens);
      set(gca,'YScale','log');
      logticks(0,1);
      set(gca,'XTick',1:length(lbl));
      set(gca,'XTickLabel',lbl);
      set(gca,'XTickLabelRotation',90);
      ylabel('Sensitivity');
      title('File Sensitivity');
      % Check sensitivity by target
      for i=1:length(obj.names)
        sel=obj.contains(i,:);
        tsens(i)=nanmedian(obj.ic(i,sel)./sens(sel));
      end
      setfig('Target Sensitivity');clf;
      bar(tsens);
      set(gca,'YScale','log');
      logticks(0,1);
      ticks=1:80:length(tsens);
      set(gca,'XTick',ticks)
      set(gca,'XTickLabel',obj.names(ticks));
      set(gca,'XTickLabelRotation',90);
      ylabel('Sensitivity (Ion count for ref file)');
      title('Target Sensitivity');

      obj.normic=obj.ic;
      for i=1:size(obj.normic,1)
        obj.normic(i,:)=obj.normic(i,:)/tsens(i);
      end
      for i=1:size(obj.normic,2)
        obj.normic(:,i)=obj.normic(:,i)/sens(i);
      end
      obj.tsens=tsens;
      obj.fsens=sens;
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
      defaults=struct('mzdata',[],'adduct',[]);
      args=processargs(defaults,varargin);

      if ischar(name)
        ind = obj.find(name,args.adduct);
      else
        ind=name;
      end
      if length(ind)>1
        % Show all adducts
        for i=1:length(ind)
          obj.getinfo(name,'adduct',obj.adduct{ind(i)});
          fprintf('\n');
        end
        return;
      end
      
      meanic=nanmean(obj.ic(ind,obj.contains(ind,:)));
      minic=nanmin(obj.normic(ind,obj.contains(ind,:)));
      meant=nanmean(obj.time(ind,obj.contains(ind,:)));
      fprintf('%s[%s] (%d): m/z=%8.4f t=%7.2f meanic=%.0f sens=%.0f\n',obj.names{ind},obj.adduct{ind}, ind, obj.mztarget(ind),meant,meanic,obj.tsens(ind));
      isomers=setdiff(find(abs(obj.mztarget-obj.mztarget(ind))<obj.MZFUZZ*2),ind);
      if length(isomers)>0
        fprintf('Isomers:\n');
        for ii=1:length(isomers)
          i=isomers(ii);
          imeanic=nanmean(obj.ic(i,obj.contains(i,:)));
          imeant=nanmean(obj.time(i,obj.contains(i,:)));
          fprintf('\t%s[%s] (%d): m/z=%8.4f (d=%.0f) t=%7.2f (d=%.0f) meanic=%.0f\n',obj.names{i},obj.adduct{i}, i, obj.mztarget(i),(obj.mztarget(i)-obj.mztarget(ind))*1e4,imeant,imeant-meant,imeanic);
        end
      end
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
        [~,filename]=fileparts(obj.files{j});
        fprintf('%-15.15s: sens=%4.2f, m/z=%8.4f t=%4.0f ic=%8.0f(%8.3f)',filename,obj.fsens(j),obj.mz(ind,j),obj.time(ind,j),obj.ic(ind,j),obj.normic(ind,j));
        m=obj.multihits{ind,j};
        if ~isempty(m)
          for k=1:length(m.mz)
            fprintf(' [mz=%8.4f, t=%4.0f, ic=%8.0f]', m.mz(k), m.time(k), m.ic(k));
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
        thresh=minic/2;
        fprintf('False positives with  m/z in [%.3f,%.3f], T in [%.0f,%.0f], NormIC >= %.3f:\n',...
                obj.mztarget(ind)+obj.MZFUZZ*[-1,1],...
                nanmean(obj.time(ind,obj.contains(ind,:)))+obj.TIMEFUZZ*[-1,1],...
                thresh);
        for j=1:length(obj.files)
          [~,filename]=fileparts(obj.files{j});
          if ~obj.contains(ind,j) && obj.normic(ind,j)>=thresh
            fprintf('%-15.15s: m/z=%8.4f t=%4.0f ic=%8.0f(%8.3f)\n',filename,obj.mz(ind,j), obj.time(ind,j), obj.ic(ind,j),obj.normic(ind,j));
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
      normic=obj.normic;
      normic(obj.contains)=0;
      [~,ord]=sort(normic(:),'desc','MissingPlacement','last');
      for i=1:args.nlist
        [ii,ij]=ind2sub(size(obj.contains),ord(i));
        [~,fname]=fileparts(obj.files{ij});
        fprintf('%12.12s:%-7.7s IC=%8.0f(%8.3f)', fname, obj.names{ii}, ic(ord(i)),normic(ord(i)));
        alias=find(obj.contains(:,ij) & abs(obj.mztarget'-obj.mztarget(ii))<=obj.MZFUZZ & abs(obj.meantime-obj.meantime(ii))<=obj.TIMEFUZZ);
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
