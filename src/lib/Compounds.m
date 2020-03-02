% Data structure to hold information about compounds located in mass spec runs
classdef Compounds < handle
  properties
    names;   % names{i} - Name of compound i (e.g.'51B07')
    mass;    % mass(i) - monoisotpic mass of compounds i
    files;   % files{j} - Mass spec filename j
    moles;   % moles(j) - Moles of each compound loaded in run j
    group;   % group{j} - name of group this file belongs to
    contains;% contains(i,j) is true if we expected to find compound i in file j
             % For compounds uniquely located in a particular run:
    mz;      % mz(i,k,j) contains the observed m/z peak on compound i, from file j for adduct k
    time;    % time(i,k,j) contains the elution time on compound i, from file j for adduct k
    meantime;% meantime(i) is the computed mean elution time of compounds i
    filetime;% filetime(i,k,j) contains the elution time on compound i, from file j without any time remapping
    ic;      % ic(i,k,j) contains the total ion count for compound i, from file j with adduct k
    normic;  % Normalized ion count (by file and by target) = ic(i,j)/tsens(i)/fsens(j)
    multihits;% multihits{i,k,j} is the list of all peaks for target i in file j with adduct k
    tsens;    % tsens(i) is the relative sensitivity to target i
    fsens;    % fsens(j) is the relative sensitivity for file j
  end
  
  properties(Constant)
    MZFUZZ=0.004; % 0.003;
    TIMEFUZZ=50;   % in seconds
    ADDUCTS=struct('name',{'M+H','M+Na','M+K'},'mass',{1.007825035,22.9897677,38.963708});
  end
  
  methods
    function obj=Compounds()
      obj.multihits={};
      obj.contains=false(0,0);
    end

    function ind=find(obj,name)
    % Find given name
    % Try both format stored in names (e.g. 'CDIV0051-B07') and short format (e.g. '51B07')
      ind=find(strcmp(name,obj.names));
      if length(ind)==0
        error('Compound %s not found',name);
      end
    end
      
    function findex=lookupMS(obj,ms)
    % Find the index of a file by name, create if missing
      if ismember(ms.path,obj.files)
        findex=find(strcmp(ms.path,obj.files));
      else
        obj.files{1,end+1}=ms.path;
        findex=length(obj.files);
        obj.mz(:,:,findex)=nan;
        obj.time(:,:,findex)=nan;
        obj.filetime(:,:,findex)=nan;
        obj.ic(:,:,findex)=nan;
        obj.normic(:,:,findex)=nan;
        obj.contains(:,findex)=false;
        obj.moles(findex)=ms.moles;
        obj.fsens(findex,1:length(obj.ADDUCTS))=nan;
      end
    end

    function [p,r,c]=getposition(obj,i)
    % Get plate (str),row (int),col(int) of compound i
      A=sscanf(obj.names{i},'%d%c%2d');
      p=A(1);
      r=A(2)-'A'+1;
      c=A(3);
    end
    
    function mz=mztarget(obj,i,k)
    % Get m/z target for compound i, adduct k
      if nargin<3
        k=[];
      end
      if nargin<2
        i=[];
      end
      if isempty(i) && isempty(k)
        mz=obj.mass'+[obj.ADDUCTS.mass];
      elseif isempty(i) 
        mz=obj.mass+obj.ADDUCTS(k).mass;
      elseif isempty(k) 
        mz=obj.mass(i)+[obj.ADDUCTS.mass];
      else
        mz=obj.mass(i)+obj.ADDUCTS(k).mass;
      end
    end
    
    function allid=checkComposition(obj,ms,varargin)  % TODO - test
    % Check composition in mass spec file using current set of compounds
      defaults=struct('debug',false,'map',[]);  
      args=processargs(defaults,varargin);
      if isempty(args.map)
        args.map=struct('mz',[0 0; 1 1 ],'time',[0 0; 1 1 ]);
      end
      allid=[];
      for k=1:length(obj.ADDUCTS)
        id=[];
        for i=1:length(obj.names)
          id=[id,struct('name',obj.names(i),'adduct',obj.ADDUCTS(k).name,'mztarget',obj.mztarget(i,k),'desc','','findargs','','mz',nan,'time',nan,'ic',nan,'relic',nan)];
          if isfinite(obj.meantime(i))
            % Only look for compounds with a known elution time
            isomers=find(abs(obj.meantime-obj.meantime(i))<obj.TIMEFUZZ & any(abs(obj.mztarget(i,k)-obj.mztarget)<obj.MZFUZZ,2));
            mztargetMS=interp1(args.map.mz(:,1),args.map.mz(:,2),obj.mztarget(i,k),'linear','extrap');
            idtmp=ms.findcompound(mztargetMS,'elutetime',obj.meantime(i),'timetol',obj.TIMEFUZZ,'mztol',obj.MZFUZZ,'debug',args.debug);
            id(end).findargs=idtmp.findargs;
            if ~isempty(idtmp.ic)
              id(end).ic=sum(idtmp.ic);
              id(end).time=mean(idtmp.time);
              id(end).mz=mean(idtmp.mz);
            end
            if length(isomers)>1 && sum(id(end).ic)>0
              fprintf('Isomer at compounds %s with ion count = %.0f\n', sprintf('%d,',isomers),sum(id(end).ic));
              id(end).ic=nan;   % Can't use it
            end
            id(end).relic=id(i).ic/obj.tsens(i,k);
          end
        end
        allid=[allid,id'];
      end
    end
    
    function [matic,id,refid]=plotComposition(obj,ms,varargin)   % TODO - fix
      defaults=struct('debug',false,'thresh',0.02,'ref',[],'adduct',1,'map',[]);  
      args=processargs(defaults,varargin);

      id=obj.checkComposition(ms,'debug',args.debug,'map',args.map);
      id=id(:,args.adduct);
      if ~isempty(args.ref)
        refid=obj.checkComposition(args.ref,'debug',args.debug);
        refid=refid(:,args.adduct);
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
      
      fprintf('%s: Located %d compounds with relative ion count >%.2f, %d with >0, out of %d with known elution time, %d total compounds\n', ms.name, sum(ic>=args.thresh), args.thresh, sum(ic>0), sum(isfinite(obj.meantime)), length(ic));
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
      
    function addCompound(obj,name,mass)
    % Add a compound
      if any(strcmp(obj.names,name))
        error('%ss already added', name);
      end

      obj.names{end+1}=name;
      nindex=length(obj.names);
      obj.mass(nindex)=mass;
      obj.tsens(nindex,1:length(obj.ADDUCTS))=nan;
      if length(obj.files)>0
        obj.mz(nindex,:,:)=nan;
        obj.time(nindex,:,:)=nan;
        obj.meantime(nindex)=nan;
        obj.filetime(nindex,:,:)=nan;
        obj.ic(nindex,:,:)=nan;
        obj.normic(nindex,:,:)=nan;
        obj.contains(nindex,:)=false;
        obj.multihits{nindex,:,:}=[];
      else
        obj.mz=nan(nindex,length(obj.ADDUCTS), 0);
        obj.time=nan(nindex,length(obj.ADDUCTS), 0);
        obj.meantime=nan(nindex,1);
        obj.filetime=nan(nindex,length(obj.ADDUCTS), 0);
        obj.ic=nan(nindex,length(obj.ADDUCTS), 0);
        obj.normic=nan(nindex,length(obj.ADDUCTS), 0);
        obj.contains=false(nindex,0);
        obj.multihits=cell(nindex,length(obj.ADDUCTS),0);
      end
    end
    
    function addCompoundsFromSDF(obj,sdf)
    % Add all the compounds in the given SDF file using a name formed from the PLATE and WELL
      for i=1:length(sdf.sdf)
        s=sdf.sdf(i);
        name=sprintf('%d%s',str2num(s.BATCH_PLATE(5:end)),s.BATCH_WELL);
        obj.addCompound(name,s.MonoisotopicMass);
      end
    end
    
    function assignTimes(obj,varargin)
    % Collect all the MS data and assign elution times to compounds
    % For each compound, 
    %   build list of observations across massspec files (elution time,whether compound is expected)
    %   find elution time with highest correlation of hit vs. expected
      defaults=struct('debug',false,'timetol',obj.TIMEFUZZ,'minhits',3,'minhitfrac',0.6);
      args=processargs(defaults,varargin);

      fprintf('assignTimes:\n');
      obj.ic(:)=nan;
      obj.normic(:)=nan;
      obj.mz(:)=nan;
      obj.time(:)=nan;
      obj.filetime(:)=nan;
      obj.meantime(:)=nan;
      for k=1:length(obj.ADDUCTS)
        % Try assignment for each adduct; only work on ones that haven't been assigned using a prior adduct
        fprintf('\n\n[%s] ',obj.ADDUCTS(k).name);
        for i=1:length(obj.names)
          if k>1 && isfinite(obj.meantime(i))
            continue;
          end
          fprintf('\n%s: ',obj.names{i});
          etimes=[];ic=[];srcadduct=[]; srcfile=[];
          for j=1:length(obj.files)
            m=obj.multihits{i,k,j};
            if ~isempty(m)
              etimes=[etimes,m.time];
              ic=[ic,m.ic];
              srcadduct=[srcadduct,repmat(k,1,length(m.time))];
              srcfile=[srcfile,repmat(j,1,length(m.time))];
            end
          end
          cont=obj.contains(i,srcfile);
          if length(etimes(cont))>1
            esort=sort(etimes(cont));
            fprintf('esort=[%s]\n',sprintf('%.0f ',esort));
            best=[1,1];bestscore=-1e10;besttime=nan;
            falseweight=length(unique(srcfile(~cont)))/length(unique(srcfile(cont)));
            for m=1:length(esort)
              for n=ceil(m+max(bestscore-2,.001)-1):length(esort)
                if esort(n)-esort(m) > 2*args.timetol
                  break
                end
                meantime=median(esort(m:n));
                % Make sure we include m,n in our tolerances
                if meantime-args.timetol > esort(m)
                  meantime=esort(m)+args.timetol;
                elseif meantime+args.timetol < esort(n)
                  meantime=esort(n)-args.timetol;
                end
                sel=abs(etimes-meantime)<=args.timetol;
                nfalse=length(unique(srcfile(~cont & sel)));
                ntrue=length(unique(srcfile(cont & sel)));
                
                %fprintf('m=%d,n=%d,true=%d,false=%d,width=%.0f,rng=[%.0f,%.0f]\n',m,n,ntrue,nfalse,esort(n)-esort(m),esort([m,n]));
                if ntrue-nfalse/falseweight>bestscore|| (ntrue-nfalse/falseweight==bestscore && esort(n)-esort(m)<esort(best(2))-esort(best(1)))
                  best=[m,n];
                  bestscore=ntrue-nfalse/falseweight;
                  besttime=meantime;
                  %fprintf('bestcore=%.1f\n',bestscore);
                end
              end
            end
            rng=besttime+args.timetol*[-1,1];
            nfalse=length(unique(srcfile(~cont & abs(etimes-besttime)<=args.timetol)));
            ntrue=length(unique(srcfile(cont & abs(etimes-besttime)<=args.timetol)));
            fprintf('Best has %d true and %d false hits over [%.0f,%.0f], width=%.0f\n',ntrue,nfalse,rng,diff(esort(best)));
          elseif length(etimes)==1
            besttime=etimes;
            ntrue=length(etimes(cont));
            nfalse=length(etimes(~cont));
          else
            fprintf('no hits\n');
            continue;
          end
          nmax=sum(obj.contains(i,:));
          fprintf('%d/%d hits with %d false hits at T=%.0f ', ntrue, nmax, nfalse, besttime);
          if ntrue/nmax < args.minhitfrac 
            fprintf('too few hits (%d hits/%d files < %.2f)',ntrue,nmax,args.minhitfrac);
          elseif ntrue<args.minhits
            fprintf('too few hits (%d hits < %d)',ntrue,args.minhits);
          else
            % Construct consensus view
            obj.meantime(i)=besttime;
            for kk=1:length(obj.ADDUCTS)
              for j=1:length(obj.files)
                m=obj.multihits{i,kk,j};
                if ~isempty(m)
                  sel=abs(m.time-besttime)<args.timetol;
                  if sum(sel)>0
                    obj.ic(i,kk,j)=sum(m.ic(sel));
                    obj.mz(i,kk,j)=sum(m.mz(sel).*m.ic(sel))/obj.ic(i,kk,j);
                    %obj.filemz(kk,j)=sum(m.filemz(sel).*m.ic(sel))/obj.ic(kk,j);
                    obj.time(i,kk,j)=sum(m.time(sel).*m.ic(sel))/obj.ic(i,kk,j);
                    obj.filetime(i,kk,j)=sum(m.filetime(sel).*m.ic(sel))/obj.ic(i,kk,j);
                  end
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
            plot((besttime-args.timetol)*[1,1],ax(3:4),':g');
            plot((besttime+args.timetol)*[1,1],ax(3:4),':g');
            xlabel('Elution time');
            ylabel('Ion count');
            title(obj.names{i});
            figure(gcf);
            keyboard;
          end
        end
        fprintf('Have times for %d compounds\n', sum(isfinite(obj.meantime)));
      end
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
          fprintf('addMS: contains list has %d nonexistent compounds names\n',length(badnames));
        end
        obj.contains(:,findex)=ismember(obj.names,args.contains);
      end

      % Attempt to locate each one uniquely
      maxic=[];   % Maximum IC per target
      for i=1:length(obj.mass)
       for k=1:length(obj.ADDUCTS)
        mztargetMS=interp1(args.map.mz(:,1),args.map.mz(:,2),obj.mztarget(i,k),'linear','extrap');
        id=ms.findcompound(mztargetMS,'mztol',args.mztol,'timetol',args.timetol,'debug',args.debug);
        % Map M/Z, time back to global values
        id.filemz=id.mz;
        id.filetime=id.time;
        id.mz=interp1(args.map.mz(:,2),args.map.mz(:,1),id.mz,'linear','extrap');
        id.time=interp1(args.map.time(:,2),args.map.time(:,1),id.time,'linear','extrap');
        obj.multihits{i,k,findex}=id;
        maxic(i,k)=max([0,id.ic]);
       end
      end
      for k=1:length(obj.ADDUCTS)
        ice=maxic(obj.contains(:,findex),k);
        icu=maxic(~obj.contains(:,findex),k);
        if isempty(icu) || nanmedian(icu(:))<200
          minic=prctile(ice(ice(:)>0),10);
        else
          minic=sqrt(median(ice(ice(:)>0))*median(icu(icu(:)>0)));
        end
        p=[25,50,75];
        fprintf('%5s: Expected have IC=[%s], unexpected have IC=[%s] (%s) @%.0f: %.1f%%,%.1f%% \n', obj.ADDUCTS(k).name, sprintf('%.0f ',prctile(ice(:),p)), sprintf('%.0f ',prctile(icu(:),p)), sprintf('%d%% ',p),minic,100*mean(ice(:)>minic),100*mean(icu(:)>minic));
        nhits=nan(length(obj.mass),1);
        for i=1:length(obj.mass)
          nhits(i)=sum(obj.multihits{i,k,findex}.ic>minic);
        end
        fprintf('       Have hits for %d/%d (with %d unique) expected compounds and %d unexpected ones with IC>=%.0f\n', sum(obj.contains(:,findex) & nhits>0), sum(obj.contains(:,findex)), sum(obj.contains(:,findex) & nhits==1), sum(~obj.contains(:,findex)&nhits>0),minic);
      end
    end
    
    function summary(obj)
    % Summarize data available
      fprintf('summary():\n');
      fprintf('Contains %d files, %d compounds, %d adducts (%d with elution time)\n', length(obj.files), length(obj.names), length(obj.ADDUCTS), sum(isfinite(obj.meantime)));
      for i=1:length(obj.files)
        [~,filename,~]=fileparts(obj.files{i});
        fprintf('%2d %-20.20s %3d/%3d/%3d compounds identified/isolated/total\n', i,filename, sum(obj.contains(:,i) & any(isfinite(obj.mz(:,:,i)),2)),sum(obj.contains(:,i)&isfinite(obj.meantime)), sum(obj.contains(:,i)));
      end
    end
    
    function notfound(obj)
    % Show compounds that have been isolated, but not found in files where they should be
      for i=1:length(obj.files)
        [~,filename,~]=fileparts(obj.files{i});
        fprintf('%2d %-20.20s %3d/%3d/%3d missing: ', i,filename, sum(obj.contains(:,i) & any(isfinite(obj.mz(:,:,i)),2)),sum(obj.contains(:,i)&isfinite(obj.meantime)), sum(obj.contains(:,i)));
        ind=obj.contains(:,i)&isfinite(obj.meantime)&~any(isfinite(obj.mz(:,:,i)),2);
        fprintf('%s\n',strjoin(obj.names(ind),','));
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
        x(ii).mass=obj.mass(i);
        
        for j=1:length(ugroups)
         files=find(strcmp(obj.group,ugroups{j})& obj.contains(i,:));
         if ~isempty(files)
         for k=1:length(obj.ADDUCTS)
          x(ii).mzoffset(j,k)=nanmean(obj.mz(i,k,files),3)-obj.mztarget(i,k);
          x(ii).elution(j,k)=nanmean(obj.time(i,k,files),3);
          x(ii).ioncount(j,k)=nanmean(obj.ic(i,k,files),3);
          f='';
          for k=1:length(files)
            [~,filename,~]=fileparts(obj.files{files(k)});
            f=[f,filename,','];
          end
          f=f(1:end-1);  % Remove trailing comma
          x(ii).files{j}=f;
          end
        end
       end
      end
      x=struct2table(x);
    end
    
    function pcolorplot(obj,adduct)
    % Heat map of ioncount in matrix of compounds * files
      if nargin<2
        adduct=1;
      end
      data=squeeze(obj.ic(:,adduct,:));
      % Normalize for xxth percentage
      for i=1:size(data,2)
        data(:,i)=data(:,i)/prctile(data(:,i),80);
      end
      data(data<.01)=0;
      data(data>10)=10;
      ti=['compounds ',obj.ADDUCTS(adduct).name];
      setfig(ti);clf;
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
      title(ti);
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
      rng=[min(min(obj.mztarget)),max(max(obj.mztarget))];
      fprintf('File     Added Offset(*1e4) [%5.0f, %5.0f]\n',rng);
      for j=1:length(udirs)
        nexttile;
        allx=[];ally=[];
        h=[];
        leg={};
        for i=1:length(obj.files)
          if strcmp(dirs{i},udirs{j})
            x=obj.mztarget;y=obj.mz(:,:,i);
            err=nanmedian(x(:)-y(:));
            fit=robustfit(x(:),y(:));
            [~,filename]=fileparts(obj.files{i});
            fprintf('%-20.20s  %5.1f [%5.1f, %5.1f]\n',filename, err*1e4,1e4*(fit(1)+(fit(2)-1)*rng));
            h(end+1)=plot(x(:),1e4*(y(:)-x(:)),'o');
            hold on;
            plot(rng,1e4*(fit(1)+(fit(2)-1)*rng),'-','Color',get(h(end),'Color'));
            leg{end+1}=filename;
            allx=[allx;x(:)];
            ally=[ally;y(:)];
          end
        end
        sel=strcmp(dirs,udirs{j});
        fit=robustfit(allx,ally);
        [~,dirname]=fileparts(udirs{j});
        fprintf('%-20.20s  %5.1f [%5.1f, %5.1f] over %s\n', '', nanmedian(allx-ally)*1e4, 1e4*(fit(1)+(fit(2)-1)*rng),dirname);
        h(end+1)=plot(rng,1e4*(fit(1)+(fit(2)-1)*rng),'k','linewidth',2);
        plot(rng,1e4*(fit(1)+(fit(2)-1)*rng+obj.MZFUZZ),'k:');
        plot(rng,1e4*(fit(1)+(fit(2)-1)*rng-obj.MZFUZZ),'k:');
        leg{end+1}='All';
        legend(h,leg,'location','best');
        xlabel('True M/Z');
        ylabel('File-True M/Z *1e4');
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
      for i=1:length(obj.names)
       for k=1:length(obj.ADDUCTS)
        for j=1:length(obj.files)
          if obj.contains(i,j)
            m=obj.multihits{i,k,j};
            if ~isempty(m) && ~isempty(m.ic)
              [ic(i,k,j),ord]=max(m.ic);
              t(i,k,j)=m.filetime(ord);
            end
          end
        end
       end
      end
      dirs={};
      for i=1:length(obj.files)
        dirs{i}=fileparts(obj.files{i});
      end
      udirs=unique(dirs);

      tref=t(:,:,ref);tref=tref(:);
      trefs=sort(tref);
      setfig('checktime');clf;
      tiledlayout('flow');
      map=[];
      for j=1:length(udirs);
        nexttile;
        h=[];leg={};
        for i=1:length(obj.files)
          if strcmp(dirs{i},udirs{j})
            t2=t(:,:,i);t2=t2(:);
            if sum(isfinite(t2))<4
              continue;
            end
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
        tmed=nanmedian(t(:,:,sel),3);
        fit=piecewise(tref,tmed(:),args.timetol,2);
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
      obj.normic=obj.ic;

      for k=1:length(obj.ADDUCTS)
      sens=[]; lbl={};
      for i=1:length(obj.files)
        sel=all(obj.contains(:,[i,ref]),2);
        sens(i)=nanmedian(obj.ic(sel,k,i)./obj.ic(sel,k,ref));
        [~,lbl{i}]=fileparts(obj.files{i});
      end
      ti=['File Sensitivity - ',obj.ADDUCTS(k).name];
      setfig(ti);clf;
      bar(sens);
      set(gca,'YScale','log');
      logticks(0,1);
      set(gca,'XTick',1:length(lbl));
      set(gca,'XTickLabel',lbl);
      set(gca,'XTickLabelRotation',90);
      ylabel('Sensitivity');
      title(ti);
      % Check sensitivity by target
      for i=1:length(obj.names)
        sel=obj.contains(i,:);
        tsens(i)=nanmedian(squeeze(obj.ic(i,k,sel))./sens(sel)');
      end
      ti=['Target Sensitivity - ',obj.ADDUCTS(k).name];
      setfig(ti);clf;
      bar(tsens);
      set(gca,'YScale','log');
      logticks(0,1);
      ticks=1:80:length(tsens);
      set(gca,'XTick',ticks)
      set(gca,'XTickLabel',obj.names(ticks));
      set(gca,'XTickLabelRotation',90);
      ylabel('Sensitivity (Ion count for ref file)');
      title(ti);

      for i=1:size(obj.normic,1)
        obj.normic(i,k,:)=obj.ic(i,k,:)/tsens(i);
      end
      for i=1:size(obj.normic,2)
        obj.normic(:,k,i)=obj.normic(:,k,i)/sens(i);
      end
      obj.tsens(:,k)=tsens;
      obj.fsens(:,k)=sens;
      end
    end
    
    function plotcompare(obj,k,f1,f2) % TODO - broken
    % Plot comparison of each compound the occurs in both f1 and f2
      [~,fname1]=fileparts(obj.files{f1});
      [~,fname2]=fileparts(obj.files{f2});
      sel=obj.contains(:,f1) & obj.contains(:,f2);
      mztarget=obj.mztarget([],k);
      mztarget=mztarget(sel);
      ti=sprintf('%s vs %s',fname1,fname2);
      setfig(ti);clf;
      subplot(221)
      plot(obj.mz(sel,k,f1)-mztarget',obj.mz(sel,k,f2)-mztarget','o');
      hold on;
      ax=axis;
      plot(ax(1:2),obj.MZFUZZ*[1,1],'r:');
      plot(ax(1:2),-obj.MZFUZZ*[1,1],'r:');
      ax=axis;
      plot(obj.MZFUZZ*[1,1],ax(3:4),'r:');
      plot(-obj.MZFUZZ*[1,1],ax(3:4),'r:');
      xlabel(fname1,'Interpreter','none');
      ylabel(fname2,'Interpreter','none');
      title('m/z offset');

      subplot(222);
      plot(mean(obj.time(sel,k,[f1,f2]),3),diff(obj.time(sel,k,[f1,f2]),[],3),'o');
      hold on;
      ax=axis;
      plot(ax(1:2),obj.TIMEFUZZ*[1,1],'r:');
      plot(ax(1:2),-obj.TIMEFUZZ*[1,1],'r:');
      xlabel('Mean (s)');
      ylabel('Diff (s)');
      title('Elution Times');
      
      ratio=obj.ic(sel,k,f2)./obj.ic(sel,k,f1);
      scaling=nanmedian(ratio(ratio>0));
      ratio=ratio/scaling;
      ratio(ratio==0)=.001;
      ratio(obj.ic(sel,k,f1)==0 & obj.ic(sel,k,f2)~=0)=1000;

      subplot(223);
      x=obj.ic(sel,k,f1); y=obj.ic(sel,k,f2);
      x(x==0)=1;
      y(y==0)=1;
      loglog(x,y,'o');
      hold on;
      ax=axis;
      plot(ax(1:2),ax(1:2)*scaling,':');
      xlabel(fname1,'Interpreter','none');
      ylabel(fname2,'Interpreter','none');
      title('Ion Counts');
      
      subplot(224);
      %      semilogy(obj.mztarget,max(ratio,.01),'o');
      histogram(log10(ratio),20);
      xlabel('log10(IC2/IC1)');
      ylabel('N');
      title(sprintf('Ion Count Ratio (scaling=%.2f)',scaling));
      
      h=suptitle(sprintf('%s N=%d',ti,sum(sel)));
      set(h,'Interpreter','none');
    end

    function plotmap(obj,k)
    % Plot map of compounds in mass vs time space
      setfig(sprintf('Compound Map [%s]',obj.ADDUCTS(k).name));clf;
      etime=[];
      for i=1:length(obj.mass)
        etime(i,1)=nanmin(obj.time(i,k,:));
        etime(i,2)=nanmax(obj.time(i,k,:));
        mz(i,1)=nanmin(obj.mz(i,k,:));
        mz(i,2)=nanmax(obj.mz(i,k,:));
      end
      plot(obj.mztarget([],k),obj.meantime,'.b');
      hold on;
      ngood=0;
      for i=1:size(etime,1)
        if isnan(obj.meantime(i))
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
      title(sprintf('Compound Map [%s] (%d distinguishable/%d identified)',obj.ADDUCTS(k).name,ngood,sum(isfinite(obj.meantime))));
    end
    
    function getinfo(obj,name,varargin)
      defaults=struct('mzdata',[],'adduct',[]);
      args=processargs(defaults,varargin);

      if ischar(name)
        ind = obj.find(name);
      else
        ind=name;
      end
      if isempty(args.adduct)
        % Show all adducts
        for i=1:length(obj.ADDUCTS)
          obj.getinfo(name,'adduct',i);
          fprintf('\n');
        end
        return;
      end
      k=args.adduct;
      
      meanic=nanmean(obj.ic(ind,k,obj.contains(ind,:)));
      minic=nanmin(obj.normic(ind,k,obj.contains(ind,:)));
      fprintf('%s[%s] (%d): m/z=%8.4f t=%7.2f sens=%.0f\n',obj.names{ind},obj.ADDUCTS(k).name, ind, obj.mztarget(ind,k),obj.meantime(ind),obj.tsens(ind,k));
      isomers=setdiff(find(abs(obj.mass-obj.mass(ind))<obj.MZFUZZ*2),ind);
      % TODO: Fix to handle same m/z (with adducts) instead of same mass
      if length(isomers)>0
        fprintf('Isomers:\n');
        for ii=1:length(isomers)
          i=isomers(ii);
          imeanic=nanmean(obj.ic(i,obj.contains(i,:)));
          fprintf('\t%s[%s] (%d): m/z=%8.4f (d=%.0f) t=%7.2f (d=%.0f) meanic=%.0f\n',obj.names{i},obj.ADDUCTS(k).name, i, obj.mass(i)+obj.ADDUCTS(k).mass,(obj.mass(i)-obj.mass(ind))*1e4,obj.meantime(i),obj.meantime(i)-obj.meantime(ind),imeanic);
        end
      end
      if ~isempty(args.mzdata)
        setfig(obj.names{ind});
        t=tiledlayout('flow');
        title(t,sprintf('%s m/z=%.4f t=%.0f',obj.names{ind},obj.mztarget(ind,k),obj.meantime(ind)));
      end
      
      for j=1:length(obj.files)
        if ~obj.contains(ind,j)
          % TODO: Could list false positives here
          continue;
        end
        [~,filename]=fileparts(obj.files{j});
        fprintf('%-15.15s: sens=%4.2f, m/z=%8.4f t=%4.0f ic=%8.0f(%8.3f)',filename,obj.fsens(j,k),obj.mz(ind,k,j),obj.time(ind,k,j),obj.ic(ind,k,j),obj.normic(ind,k,j)/obj.fsens(j,k));
        m=obj.multihits{ind,k,j};
        if ~isempty(m)
          for p=1:length(m.mz)
            if isnan(obj.time(ind,k,j)) || abs(m.time(p)-obj.time(ind,k,j))>1
              fprintf(' [mz=%8.4f, t=%4.0f, ic=%8.0f]', m.mz(p), m.time(p), m.ic(p));
            end
          end
        end
        if ~isempty(args.mzdata)
          nexttile;
          obj.plotscan(ind,args.mzdata{j});
        end
        fprintf('\n');
      end
      if isfinite(meanic) && isfinite(obj.meantime(ind))
        % False positives
        thresh=minic/2;
        fprintf('False positives with  m/z in [%.3f,%.3f], T in [%.0f,%.0f], NormIC >= %.3f:\n',...
                obj.mztarget(ind,k)+obj.MZFUZZ*[-1,1],...
                obj.meantime(ind)+obj.TIMEFUZZ*[-1,1],...
                thresh);
        for j=1:length(obj.files)
          [~,filename]=fileparts(obj.files{j});
          if ~obj.contains(ind,j) && obj.normic(ind,k,j)>=thresh
            fprintf(' %-14.14s: sens=%4.2f, m/z=%8.4f t=%4.0f ic=%8.0f(%8.3f)\n',filename,obj.fsens(j,k),obj.mz(ind,k,j), obj.time(ind,k,j), obj.ic(ind,k,j),obj.normic(ind,k,j));
            if ~isempty(args.mzdata)
              nexttile;
              obj.plotscan(ind,args.mzdata{j});
            end
          end
        end
      end
    end

    function plotscan(obj,ind,mzdata)
    % TODO: Broken
      meanic=nanmean(obj.ic(ind,obj.contains(ind,:)));
      [ic,mz,t]=mzdata.mzscan(obj.mztarget(ind,k),'mztol',obj.MZFUZZ);
      plot(t,ic);
      ax=axis;
      hold on;
      plot(obj.meantime(ind)+obj.TIMEFUZZ*[1,1],ax(3:4),':b');
      plot(obj.meantime(ind)-obj.TIMEFUZZ*[1,1],ax(3:4),':b');
      ylabel('Ion Count');
      yyaxis right
      plot(t,mz,'r');
      hold on;
      ax=axis;
      axis(ax);
      plot(ax(1:2),obj.mztarget(ind,k)+obj.MZFUZZ*[1,1],':r');
      plot(ax(1:2),obj.mztarget(ind,k)-obj.MZFUZZ*[1,1],':r');
      ylabel('M/Z');
      if isfinite(obj.meantime(ind))
        ax(1:2)=obj.meantime(ind)+obj.TIMEFUZZ*2*[-1,1];
      end
      ax(3:4)=obj.mztarget(ind,k)+obj.MZFUZZ*2*[-1,1];
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
        % TODO: Fix to handle adducts
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
  