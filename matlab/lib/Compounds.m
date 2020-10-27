% Data structure to hold information about compounds located in mass spec runs
classdef Compounds < handle
  properties
    names;   % names{i} - Name of compound i (e.g.'51B07')
    mass;    % mass(i) - monoisotpic mass of compounds i
    samples;   % Name of sample (mass spec file)
    files;   % files{j} - Mass spec filename j
    maps;    % Piecewise linear maps for converting between references (col 1) and file values (col2)
    moles;   % moles(j) - Moles of each compound loaded in run j
    group;   % group{j} - name of group this file belongs to
    contains;% contains(i,j) is true if we expected to find compound i in file j
             % For compounds uniquely located in a particular run:
    mz;     % mz(i,k,j) contains the observed m/z peak on compound i, from file j for adduct k
    time;    % time(i,k,j) contains the elution time on compound i, from file j for adduct k
    meantime;% meantime(i) is the computed mean elution time of compounds i
    timewindow;  % timewindow(i,2) is the time window for compound i elution
    ic;      % ic(i,k,j) contains the total ion count for compound i, from file j with adduct k
    normic;  % normic(i,k,j) normalized ion count = ic(i,j,k)/tsens(i)/fsens(k)
    multihits;% multihits(i,k,j) is the list of all peaks for target i in file j with adduct k (index into reffeatures)
    allfeatures;% allfeatures(j) is the list of all features in file j (copied from ms feature list)
    reffeatures;% reffeatures(j) is the list of all features in file j in reference m/z,time space using map
    featureindex;  % featureindex(i,k,j) is the index of feature finally chosen
    tsens;    % tsens(i,k) is the relative sensitivity to target i with adduct k
    fsens;    % fsens(j) is the relative sensitivity for file j
    astats;   % astats(i) - struct showing setup for assigning elution time to compound i
    sdf;      % SDF data
    MZFUZZ;   % global mztol
    TIMEFUZZ; % global timetol
    ADDUCTS;  % adducts to use
  end
  
  properties(Constant)
  end
  
  methods
    function obj=Compounds(mzfuzz,timefuzz)
      obj.multihits={};
      obj.contains=false(0,0);
      obj.samples={};
      obj.allfeatures=FeatureList.empty;
      obj.reffeatures=FeatureList.empty;
      obj.astats=struct('run',{},'args',{},'adduct',{},'sel',{},'hitgood',{},'hitlow',{},'hithigh',{},'missstrong',{},'missweak',{},'FP',{},'FN',{});
      obj.sdf=SDF();
      if nargin<1
        obj.MZFUZZ=0.006;
      else
        obj.MZFUZZ=mzfuzz;
      end
      if nargin<2
        obj.TIMEFUZZ=40/60;   % in minutes
      else
        obj.TIMEFUZZ=timefuzz;
      end
      obj.ADDUCTS=struct('name',{'M+H','M+Na','M+K'},'mass',{1.007276,22.989218,38.963158});
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
        obj.ic(:,:,findex)=nan;
        obj.normic(:,:,findex)=nan;
        obj.contains(:,findex)=false;
        obj.moles(findex)=ms.moles;
        obj.fsens(findex)=nan;
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
      defaults=struct('debug',false,'map',[],'timetol',obj.TIMEFUZZ,'mztol',obj.MZFUZZ,'remap',false);  
      args=processargs(defaults,varargin);
      if args.remap
        args.map=obj.computeMap(ms);
      elseif isempty(args.map)
        fprintf('Warning: Using default map');
        args.map=struct('mz',[0 0; 1 1 ],'time',[0 0; 1 1 ]);
      end
      allid=[];
      for k=1:length(obj.ADDUCTS)
        id=[];
        for i=1:length(obj.names)
          id=[id,struct('name',obj.names(i),'adduct',obj.ADDUCTS(k).name,'mztarget',obj.mztarget(i,k),'desc','','findargs','','mz',nan,'time',nan,'ic',nan,'relic',nan,'filemz',nan,'filetime',nan)];
          if isfinite(obj.meantime(i))
            % Only look for compounds with a known elution time
            aliases=find(abs(obj.meantime-obj.meantime(i))<args.timetol & any(abs(obj.mztarget(i,k)-obj.mztarget)<args.mztol,2));
            mztargetMS=interp1(args.map.mz(:,1),args.map.mz(:,2),obj.mztarget(i,k),'linear','extrap');
            timetarget=interp1(args.map.time(:,1),args.map.time(:,2),obj.meantime(i),'linear','extrap');
            idtmp=ms.findcompound(mztargetMS,'elutetime',timetarget,'timetol',args.timetol,'mztol',args.mztol,'debug',args.debug,'peakratio',0);
            id(end).findargs=idtmp.findargs;
            if ~isempty(idtmp.ic)
              id(end).ic=sum(idtmp.ic);
              id(end).filetime=mean(idtmp.time);
              id(end).filemz=mean(idtmp.mz);
              id(end).mz=interp1(args.map.mz(:,2),args.map.mz(:,1),id(end).filemz,'linear','extrap');
              id(end).time=interp1(args.map.time(:,2),args.map.time(:,1),id(end).filetime,'linear','extrap');
              if isfield(idtmp,'fwhh')
                id(end).fwhh=interp1(args.map.time(:,2),args.map.time(:,1),[min(idtmp.fwhh(:,1)),max(idtmp(fwhh(:,2)))]);
              end
            end
            if length(aliases)>1 && sum(id(end).ic)>0
              if args.debug
                fprintf('Alias at compounds %s with ion count = %.0f\n', sprintf('%d,',aliases),sum(id(end).ic));
              end
              id(end).ic=nan;   % Can't use it
            end
            id(end).relic=id(i).ic/obj.tsens(i,k);
          end
        end
        allid=[allid,id'];
      end
    end
    
    function labelFeatures(obj,fl,varargin)
    % Attempt to label the features in fl
    % Assume this featurelist is already aligned (mapped) to compound's times and true m/z
    % Use FeatureList.maptoref() if not
      defaults=struct('debug',false,'timetol',obj.TIMEFUZZ,'mztol',obj.MZFUZZ);
      args=processargs(defaults,varargin);

      for j=1:length(obj.ADDUCTS)
        nlabel=0;
        for i=1:length(obj.names)
          if ~isfinite(obj.meantime(i))
            % Not isolated
            continue;
          end
          cname=sprintf('%s[%s]',obj.names{i},obj.ADDUCTS(j).name);
          f=fl.getbymz(obj.mass(i)+obj.ADDUCTS(j).mass,'timerange',obj.meantime(i)+[-1,1]*args.timetol,'mztol',args.mztol);
          if length(f.features)>1 
            fprintf('Have %d possible matches for %s (m/z=%.4f, t=%.2f)\n', length(f.features), cname,obj.mass(i)+obj.ADDUCTS(j).mass, obj.meantime(i));
            tmin=min([f.features.time]);
            tmax=max([f.features.time]);
            if tmax-tmin < args.timetol
              fprintf('\tAll in time range [%.2f,%.2f] < mztol; assuming largest is this peak\n',tmin,tmax);
              f.features=f.features([f.features.intensity]==max([f.features.intensity]));
            else
              for fi=1:length(f.features)
                fprintf('\tm/z=%.4f, t=%.2f, ic=%.0f\n', f.features(fi).mz, f.features(fi).time, f.features(fi).intensity);
              end
            end
          end

          if length(f.features)==1
            if ~ismember(cname,f.features(1).labels)
              f.features(1).labels{end+1}=cname;
            end
            nlabel=nlabel+1;
          end
        end
        fprintf('Added %d labels using %s adduct\n',nlabel,obj.ADDUCTS(j).name);
      end
    end
    
    function [matic,id,refid]=plotComposition(obj,ms,varargin)   % TODO - fix
      defaults=struct('debug',false,'thresh',0.02,'ref',[],'adduct',1,'map',[],'timetol',obj.TIMEFUZZ,'mztol',obj.MZFUZZ,'refmap',[]);  
      args=processargs(defaults,varargin);

      id=obj.checkComposition(ms,'debug',args.debug,'map',args.map,'timetol',args.timetol,'mztol',args.mztol);
      id=id(:,args.adduct);
      if ~isempty(args.ref)
        refid=obj.checkComposition(args.ref,'debug',args.debug,'map',args.refmap,'timetol',args.timetol,'mztol',args.mztol);
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
        plot(ax(1:2),-args.timetol*[1,1],':r');
        plot(ax(1:2),args.timetol*[1,1],':r');
        xlabel(args.ref.name,'Interpreter','none');
        ylabel('Time diff (s)','Interpreter','none');
        legend(h,{sprintf('RelIC<=%1.g',args.thresh),sprintf('RelIC>%.1g',args.thresh)},'location','best');
        title('Elution Time Compare');
        
        fprintf('IC/Ref IC:  [1,5,50,95,99-pct]: %s\n',sprintf('%.2f ',prctile(relic(sel),[1,5,50,95,99])));
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
      %obj.sdf.sdf(end+1)=struct();
      nindex=length(obj.names);
      obj.mass(nindex)=mass;
      obj.tsens(nindex,1:length(obj.ADDUCTS))=nan;
      if length(obj.files)>0
        obj.mz(nindex,:,:)=nan;
        obj.time(nindex,:,:)=nan;
        obj.meantime(nindex)=nan;
        obj.ic(nindex,:,:)=nan;
        obj.normic(nindex,:,:)=nan;
        obj.contains(nindex,:)=false;
      else
        obj.mz=nan(nindex,length(obj.ADDUCTS), 0);
        obj.time=nan(nindex,length(obj.ADDUCTS), 0);
        obj.meantime=nan(nindex,1);
        obj.ic=nan(nindex,length(obj.ADDUCTS), 0);
        obj.normic=nan(nindex,length(obj.ADDUCTS), 0);
        obj.contains=false(nindex,0);
      end
    end
    
    function addCompoundsFromSDF(obj,sdf)
    % Add all the compounds in the given SDF file using a name formed from the PLATE and WELL
      assert(isempty(obj.names));   % Otherwise, can't set obj.sdf to correspond
      obj.sdf=sdf;
      for i=1:length(sdf.sdf)
        s=sdf.sdf(i);
        name=sprintf('%d%s',str2num(s.BATCH_PLATE(5:end)),s.BATCH_WELL);
        obj.addCompound(name,s.MostAbundantMass);  % May be different from monoisotopic mass
        index=strcmp(obj.names,name);
      end
    end
    
    function map=computeMap(obj,ms,varargin)
    % Compute time mapping for given ms file
      defaults=struct('initmap',[]);
      args=processargs(defaults,varargin);

      filetime=str2num(ms.mzxml.mzXML.msRun.endTime(3:end-1));
      objtime=3300/60;   % Baseline time
      if isempty(args.initmap)
        fprintf('Normalizing file end time of %.1f to baseline of %.1f\n', filetime, objtime);
        map=struct('mz',[0,0;1,1], 'time',[0,0;objtime,filetime]);
      else
        map=args.initmap;
      end
      ti=['computeMap-',ms.name];
      setfig(ti);clf;
      niter=3;
      for iter=1:niter
        allid=obj.checkComposition(ms,'map',map,'timetol',500/60/iter,'mztol',.02/iter);
        t=[];mz=[];  % Build lists of [reference,file]
        for i=1:length(allid(:))
          if isfinite(allid(i).time)
            % First coords are in mapped values, second in uncorrected values
            t(end+1,1:2)=[allid(i).findargs.elutetime,allid(i).time];
            mz(end+1,1:2)=[allid(i).mztarget,allid(i).mz];
          end
        end
        % Undo mapping of reference t(:,1),mz(:,1)
        t(:,1)=interp1(map.time(:,2),map.time(:,1),t(:,1),'linear','extrap');
        mz(:,1)=interp1(map.mz(:,2),map.mz(:,1),mz(:,1),'linear','extrap');
        
        subplot(niter,2,iter*2-1);
        plot(t(:,1),t(:,2),'.');
        xlabel('Reference Time');  
        ylabel([ms.name,' Time - Ref'],'Interpreter','none');
        tfit=robustfit(t(:,1),t(:,2),'cauchy',.1);
        %[tfit(2),tfit(1)]=houghregress(t(:,1),t(:,2)-t(:,1));
        hold on;
        rng=[min(t(:,1)),max(t(:,1))];
        plot(rng,rng*tfit(2)+tfit(1),'r-');
        plot(rng,rng*tfit(2)+tfit(1)+obj.TIMEFUZZ,'r:');
        plot(rng,rng*tfit(2)+tfit(1)-obj.TIMEFUZZ,'r:');
        %tfit(2)=tfit(2)+1;
        title(sprintf('Iteration %d',iter));
        
        subplot(niter,2,iter*2);
        plot(mz(:,1),mz(:,2)-mz(:,1),'.');
        xlabel('Reference m/z');  
        ylabel([ms.name,' - ref m/z'],'Interpreter','none');
        title(sprintf('Iteration %d',iter));
        %[mzfit(2),mzfit(1)]=houghregress(mz(:,1),mz(:,2)-mz(:,1));
        mzfit=robustfit(mz(:,1),mz(:,2)-mz(:,1));
        hold on;
        rng=[min(mz(:,1)),max(mz(:,1))];
        plot(rng,rng*mzfit(2)+mzfit(1),'r-');
        plot(rng,rng*mzfit(2)+mzfit(1)+obj.MZFUZZ,'r:');
        plot(rng,rng*mzfit(2)+mzfit(1)-obj.MZFUZZ,'r:');
        mzfit(2)=mzfit(2)+1;
        
        % Update mapping
        map.time(:,2)=map.time(:,1)*tfit(2)+tfit(1);
        map.mz(:,2)=map.mz(:,1)*mzfit(2)+mzfit(1);
        fprintf('After iteration %d, time fit=[%f,%f], m/z fit=[%f,%f]\n',iter,tfit,mzfit);
        pause(0.1);
      end
      h=suptitle(ti);set(h,'Interpreter','none');
    end
    
    function assignTimes(obj,varargin)
    % Collect all the MS data and assign elution times to compounds
    % For each compound, 
    %   build list of observations across features (elution time,whether compound is expected) that are not already assigned
    %   find elution times which have <= given number of false negatives and falsepositives
    %   if unique, assign to compound;  if multiples display message
      defaults=struct('debug',0,'timetol',obj.TIMEFUZZ,'minhits',3,'mztol',obj.MZFUZZ,'plot','','minic',1000,'trace',[],'maxFN',0,'maxFP',0,'clear',true,'detectionThreshold',2000,'normicrange',[0.4,2.5],'usefiles',[]);
      args=processargs(defaults,varargin);

      if args.clear || isempty(obj.astats)
        arun=1;
      else
        arun=max([obj.astats.run])+1;
      end
      
      fprintf('assignTimes:\n');
      if args.clear
        obj.ic(:)=nan;
        obj.normic(:)=nan;
        obj.mz(:)=nan;
        obj.time(:)=nan;
        obj.meantime(:)=nan;
        obj.timewindow(:,:)=nan;
        obj.featureindex=nan(length(obj.names),length(obj.ADDUCTS),length(obj.samples));
        obj.astats=obj.astats([]);
      end
      for k=1:length(obj.ADDUCTS)
        % Try assignment for each adduct; only work on ones that haven't been assigned using a prior adduct
        fprintf('[%s]\n',obj.ADDUCTS(k).name);
        for i=1:length(obj.names)
          if isfinite(obj.meantime(i))
            continue;
          end
          etimes=[];intensity=[];srcfile=[];
          mztarget=obj.mztarget(i,k);
          for j=1:length(obj.files)
            if ~isempty(args.usefiles) && ~ismember(j,args.usefiles)
              continue;
            end
            mhits=obj.multihits{i,k,j};
            if ~isempty(mhits)
              used=obj.featureindex(:,:,j);
              used=used(isfinite(used(:)));
              featinds=mhits(~ismember(mhits,used));   % This is faster than using setdiff()
              if ~isempty(featinds)
                feat=obj.reffeatures(j).features(featinds);
                mzsel=abs([feat.mz]-mztarget)<=args.mztol & [feat.intensity]>=args.minic & [feat.isotope]<=1;
                if any(mzsel)
                  etimes=[etimes,feat(mzsel).time];
                  intensity=[intensity,feat(mzsel).intensity];
                  srcfile=[srcfile,repmat(j,1,length(feat(mzsel)))];
                end
              end
            end
          end
          cont=obj.contains(i,srcfile);

          bestset=[];
          nbest=0;
          if length(etimes(cont))>0
            [esort,ord]=sort(etimes(cont));
            if args.debug || strcmp(args.plot,obj.names{i})  || ismember(i,args.trace)
              fprintf('%s esort=[%s]\n',obj.names{i}, sprintf('%.2f ',esort));
            end
            for m=1:length(esort)
              for n=m+args.minhits-1:length(esort)
                if esort(n)-esort(m) > 2*args.timetol
                  break
                end
                sel=find(etimes>=esort(m) & etimes<=esort(n));
                % Keep only one peak per srcfile; the largest one
                keep=false(size(sel));
                for is=1:length(sel)
                  keep(is)=~any(srcfile(sel)==srcfile(sel(is)) & intensity(sel)>intensity(sel(is)));
                end
                sel=sel(keep);
                selexpected=sel(cont(sel));
                if length(selexpected)<args.minhits
                  continue;
                end
                expected=obj.contains(i,:);   % Files expected to contain this target
                if ~isempty(args.usefiles)
                  % Remove the unused files from expected
                  expected(setdiff(1:length(expected),args.usefiles))=false;
                end
                
                % Estimate target sensitity given this set
                ratio=intensity(selexpected)./obj.fsens(srcfile(selexpected));
                ratio=ratio(isfinite(ratio));   % Using nanmedian() is much slower
                tsens=median(ratio);
                expectedIC=zeros(size(obj.fsens));
                expectedIC(expected)=tsens*obj.fsens(expected);
                normic=zeros(size(obj.fsens));
                normic(srcfile(sel))=intensity(sel)./expectedIC(srcfile(sel));
                % Categorize each file as
                %  trueneg - unexpected and not present
                %  falsepos - unexpected but present
                %  hitgood - expected and present and in correct normic range
                %  hitlow - expected and present but low normic
                %  hithigh - expected and present but high normic
                %  missweak - expected (weakly) and not present
                %  missstrong - expected (strongly) and not present
                present=ismember(1:length(obj.samples),srcfile(sel));
                trueneg=~present&~expected;
                falsepos=present&~expected;
                lowic=normic<args.normicrange(1);
                highic=normic>args.normicrange(2);
                hit = present & expected;
                hitlow=hit & lowic; 
                hithigh=hit &highic;
                hitgood=hit & ~hitlow & ~hithigh; 
                miss=~present & expected;
                missstrong=miss&~(expectedIC<args.detectionThreshold);   % Include nan expectedIC here
                missweak=miss&expectedIC<args.detectionThreshold;

                assert(all(trueneg|falsepos|hitgood|hitlow|hithigh|missweak|missstrong));
                assert(sum(trueneg)+sum(falsepos)+sum(hitgood)+sum(hitlow)+sum(hithigh)+sum(missweak)+sum(missstrong)==length(obj.samples));
                nFP=sum(falsepos);
                nFN=sum(missstrong|hitlow|hithigh);
                if ismember(i,args.trace)
                  fprintf('m=%d,n=%d, keep=%s, hits=(%d good,%d low, %d high), miss=(%d strong,%d weak), FP=%d,FN=%d, width=%.2f,rng=[%.2f,%.2f]\n',m,n,sprintf('%d ',sel),sum(hitgood),sum(hitlow),sum(hithigh),sum(missstrong),sum(missweak),nFP,nFN,esort(n)-esort(m),esort([m,n]));
                  if m==3 && n==11
                    keyboard;
                  end
                end
                if sum(hitgood)<args.minhits || nFP>args.maxFP || nFN > args.maxFN
                  continue;
                end
                if length(sel)>length(bestset)
                  bestset=sel;
                  astats=struct('run',arun,'args',args,'adduct',k,'sel',srcfile(sel),'hitgood',sum(hitgood),'hitlow',sum(hitlow),'hithigh',sum(hithigh),'missstrong',sum(missstrong),'missweak',sum(missweak),'FP',nFP,'FN',nFN);
                  nbest=1;
                  if ismember(i,args.trace)
                    fprintf('best set: %s\n', sprintf('%d ',sel));
                  end
                elseif length(bestset)==length(sel) && any(~ismember(sel,bestset))
                  % Different set
                  if ismember(i,args.trace)
                    fprintf('disjoint set: %s\n', sprintf('%d ',sel));
                  end
                  nbest=nbest+1;
                else
                  if ismember(i,args.trace)
                    fprintf('ignored set: %s\n', sprintf('%d ',sel));
                  end
                end
              end
            end
          end
          if isempty(bestset)
            if args.debug || ismember(i,args.trace)
              fprintf('No hits for %s with FN<=%d and FP<=%d\n', obj.names{i}, args.maxFN, args.maxFP);
            end
            continue;
          elseif nbest>1
            if args.debug  || ismember(i,args.trace)
              fprintf('Have %d equivalent elution times for %s\n', nbest, obj.names{i});
            end   
            continue;
          end

          bestwindow=[min(etimes(bestset)),max(etimes(bestset))];
          nFP=length(unique(srcfile(~cont & etimes>=bestwindow(1) & etimes<=bestwindow(2))));
          ntrue=length(unique(srcfile(cont & etimes>=bestwindow(1) & etimes<=bestwindow(2))));
          nFN=sum(obj.contains(i,:))-ntrue;
          if args.debug || strcmp(args.plot,obj.names{i})  || ismember(i,args.trace)
            fprintf('Best has %d true, %d FN and %d FP over [%.2f,%.2f]\n',ntrue,nFN,nFP,bestwindow);
          end

          % Construct consensus view
          obj.meantime(i)=mean(bestwindow);
          obj.timewindow(i,1:2)=bestwindow;
          obj.astats(i)=astats;
          for kk=1:length(obj.ADDUCTS)
            for j=1:length(obj.files)
              fi=obj.multihits{i,kk,j};
              time=[obj.reffeatures(j).features(fi).time];
              sel=find(time>=obj.timewindow(i,1) & time <= obj.timewindow(i,2));
              fi=fi(sel);
              if ~isempty(fi)
                if length(fi)>1
                  % More than one feature in time window; use closest to mean time (perhaps should use peak height?)
                  if args.debug || ismember(i,args.trace)
                    fprintf('Have %d ambiguous features for %s[%s] in %s\n', length(fi), obj.names{i},obj.ADDUCTS(k).name,obj.samples{j});
                  end
                  [~,mind]=min(abs([obj.reffeatures(j).features(fi).time]-obj.meantime(i)));
                  fi=fi(mind);
                end
                obj.featureindex(i,kk,j)=fi;
                feat=obj.reffeatures(j).features(fi);
                obj.ic(i,kk,j)=feat.intensity;
                obj.mz(i,kk,j)=feat.mz;
                obj.time(i,kk,j)=feat.time;
              end
            end
          end
          if args.debug || strcmp(args.plot,obj.names{i})  || ismember(i,args.trace)
            fprintf('\n');
          end
        
          if strcmp(args.plot,obj.names{i})
            setfig(['assignTimes-',obj.names{i}]);clf;
            semilogy(etimes(~cont),area(~cont),'or');
            hold on;
            semilogy(etimes(cont),area(cont),'xb');
            ax=axis;
            plot(bestwindow(1)*[1,1],ax(3:4),':m');
            plot(bestwindow(2)*[1,1],ax(3:4),':m');
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
    
    function remapfeatures(obj)
    % Create remapped version of file features into reference space
      fprintf('Remapping %d feature lists...',length(obj.allfeatures));
      for i=1:length(obj.allfeatures)
        if length(obj.reffeatures)<i || isempty(obj.reffeatures(i)) || length(obj.reffeatures(i).features)~=length(obj.allfeatures(i).features)
          fprintf('%d...',i);
          obj.reffeatures(i)=obj.allfeatures(i).maptoref(obj.maps{i});
        end
      end
      fprintf('done\n');
    end
    
    function addMS(obj,ms,varargin)
    % Add the unique peaks for compounds in given SDF from a particular M/S run
    % Use prior analyses to figure out the expected elution time for each compound
    % or scan all elution times if the no prior data (keep only if a unique peak is determined)
      defaults=struct('group','','contains',{{}},'map',struct('mz',[0 0; 1 1 ],'time',[0 0; 1 1 ]),'sample',[]);
      % mzmap(i,2) - piecewise linear map for M/Z; mzmap(:,1) is true M/Z, mzmap(:,2) is for values in ms file
      % timemap(i,2) - piecewise linear map for elution times; timemap(:,1) is "standard" elution time
      args=processargs(defaults,varargin);
      
      fprintf('Adding data from %s...',ms.name);
      findex=obj.lookupMS(ms);
      if ~isempty(args.sample)
        obj.samples{findex}=args.sample;
      else
        [~,filename,~]=fileparts(obj.files{findex});
        obj.samples{findex}=filename;
      end
      obj.maps{findex}=args.map;
      if ~isempty(args.group)
        obj.group{findex}=args.group;
      end
      if isempty(args.contains)
        obj.contains(:,findex)=true;
      elseif islogical(args.contains)
        assert(length(args.contains)==size(obj.contains,1));
        obj.contains(:,findex)=args.contains;
      elseif strcmp(args.contains{1},'NONE')
        obj.contains(:,findex)=false;
      else
        badnames=setdiff(args.contains,obj.names);
        if ~isempty(badnames)
          fprintf('addMS: contains list has %d nonexistent compound names\n',length(badnames));
        end
        obj.contains(:,findex)=ismember(obj.names,args.contains);
      end
      fl=ms.featurelists(end);
      obj.allfeatures(findex)=fl;
      if length(obj.reffeatures)>=findex
        % Clear it
        obj.reffeatures(findex)=FeatureList.empty;
      end
      fprintf('done\n');
    end
    
    function findfeatures(obj,varargin)
      defaults=struct('mztol',obj.MZFUZZ,'sample',1:length(obj.samples));
      args=processargs(defaults,varargin);

      obj.remapfeatures();
      for j=1:length(args.sample)
        findex=args.sample(j);
        map=obj.maps{findex};
        
        % Attempt to locate each one uniquely
        maxic=[];   % Maximum IC per target
        fprintf('%s...',obj.samples{findex});
        for i=1:length(obj.mass)
          if mod(i,100)==1
            fprintf('%d...',i);
          end
          for k=1:length(obj.ADDUCTS)
            % Use features
            [flmz,findices]=obj.reffeatures(findex).getbymz(obj.mztarget(i,k),'mztol',args.mztol);
            obj.multihits{i,k,findex}=findices;
            maxic(i,k)=max([0,flmz.features.intensity]);
          end
        end
        fprintf('done.\n');
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
            nhits(i)=sum([obj.reffeatures(findex).features(obj.multihits{i,k,findex}).intensity]>minic);
          end
          fprintf('       Have hits for %d/%d (with %d unique) expected compounds and %d unexpected ones with IC>=%.0f\n', sum(obj.contains(:,findex) & nhits>0), sum(obj.contains(:,findex)), sum(obj.contains(:,findex) & nhits==1), sum(~obj.contains(:,findex)&nhits>0),minic);
          nfeatures=cellfun(@(z) length(z),obj.multihits(:,k,findex));
          contains=obj.contains(:,findex);
          fprintf('       Have features for %d/%d=%.0f%% (with %d unique) expected compounds and %d=%.0f%% unexpected ones\n', ...
                  sum(contains & nfeatures>0), sum(contains), sum(contains&nfeatures>0)/sum(contains)*100,...
                  sum(contains & nfeatures==1), ...
                  sum(~contains&nfeatures>0),sum(~contains&nfeatures>0)/sum(~contains)*100);
        end
      end
    end
    
    function summary(obj,varargin)
      defaults=struct('thresh',0.1);
      args=processargs(defaults,varargin);
    % Summarize data available
      fprintf('summary():\n');
      fprintf('Contains %d files, %d compounds, %d adducts (%d with elution time) Counts for norm IC>=%.2f\n', length(obj.files), length(obj.names), length(obj.ADDUCTS), sum(isfinite(obj.meantime)),args.thresh);
      for i=1:length(obj.files)
        hit=any(obj.normic(:,:,i)>args.thresh,2);
        fprintf('%2d %-25.25s %3d/%3d/%3d compounds identified/isolated/total, %d false positives\n', i,obj.samples{i}, sum(obj.contains(:,i) & hit),sum(obj.contains(:,i)&isfinite(obj.meantime)), sum(obj.contains(:,i)),sum(~obj.contains(:,i) & hit));
      end
    end
    
    function notfound(obj,varargin)
      defaults=struct('minintensity',1000);   % Only show as missing if the expected intensity is more than this
      args=processargs(defaults,varargin);
    % Show compounds that have been isolated, but not found in files where they should be
      for i=1:length(obj.files)
        fprintf('%2d %-20.20s %3d/%3d/%3d missing: ', i,obj.samples{i}, sum(obj.contains(:,i) & any(isfinite(obj.mz(:,:,i)),2)),sum(obj.contains(:,i)&isfinite(obj.meantime)), sum(obj.contains(:,i)));
        ind=find(obj.contains(:,i)&isfinite(obj.meantime)&~any(isfinite(obj.mz(:,:,i)),2));
        for jj=1:length(ind)
          j=ind(jj);
          expected=max(obj.tsens(j,:))*obj.fsens(i);
          if expected>=args.minintensity
            fprintf('%s(%.0f) ',obj.names{j},expected);
          end
        end
        fprintf('\n');
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
            f=[f,obj.samples{files(k)},','];
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
      data(data>10)=10;
      ti=['compounds ',obj.ADDUCTS(adduct).name];
      setfig(ti);clf;
      t=tiledlayout(2,1);
      t.Title.String=ti;
      h=[];
      for i=1:2
        nexttile();
        pdata=data;
        if i==1
          %pdata(isnan(pdata))=0;
          pdata(~obj.contains)=nan;
          pdata(pdata<.01)=0.01;
        else
          pdata(obj.contains)=nan;
          pdata(pdata<.01)=0.01;
        end
        pdata(end+1,:)=nan;
        pdata(:,end+1)=nan;
        pcolor(log10(pdata)');
        h(i)=gca;
        shading flat;
        colorbar;
        xlabel('Compound');
        ylabel('File');
        set(gca,'YTick',(1:length(obj.files))+0.5);
        set(gca,'YTickLabels',obj.samples);
        set(gca,'ticklabelinterpreter','none');
        set(gca,'TickDir','out');
        set(gca,'TickLength',[1,1]*.002)
        ticks=1:80:length(obj.names);
        set(gca,'XTick',ticks)
        set(gca,'XTickLabel',obj.names(ticks));
        set(gca,'XTickLabelRotation',90);
        if i==1
          title('Expected');
        else
          title('Unexpected');
        end
      end
      linkaxes(h);
    end
    
    function checkmzoffset(obj,varargin)
    % Check whether mzoffset used when reading mass spec files should be changed
      defaults=struct('files',1:length(obj.files));
      args=processargs(defaults,varargin);

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
        for ii=1:length(args.files)
          i=args.files(ii);
          if strcmp(dirs{i},udirs{j})
            x=obj.mztarget;y=obj.mz(:,:,i);
            x(~obj.contains(:,i),:)=nan;   % Only ones that it is supposed to contain
            npts=sum(isfinite(y(:)) & isfinite(x(:)));
            if npts<10
              fprintf('Only %d data points for %s ... skipping\n', npts, obj.files{i});
              continue;
            end
            err=nanmedian(y(:)-x(:));
            fit=robustfit(x(:),y(:));
            fprintf('%-20.20s  %5.1f [%5.1f, %5.1f] N=%d\n',obj.samples{i}, err*1e4,1e4*(fit(1)+(fit(2)-1)*rng),npts);
            h(end+1)=plot(x(:),1e4*(y(:)-x(:)),'o');
            hold on;
            plot(rng,1e4*(fit(1)+(fit(2)-1)*rng),'-','Color',get(h(end),'Color'));
            leg{end+1}=obj.samples{i};
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
      defaults=struct('debug',false,'timetol',obj.TIMEFUZZ,'refrange',[500,2500]/60,'files',1:length(obj.files));
      args=processargs(defaults,varargin);

      args.files=union(args.files,ref);
      fprintf('checktime(%d,''timetol'',%.2f)\n',ref,args.timetol);
      t = nan(size(obj.mz));
      ic = nan(size(obj.mz));
      for i=1:length(obj.names)
       for k=1:length(obj.ADDUCTS)
        for j=args.files
          if obj.contains(i,j)
            t(i,k,j)=obj.time(i,k,j);
            ic(i,k,j)=obj.ic(i,k,j);
          end
        end
       end
      end
      dirs={};
      for i=args.files
        dirs{i}=fileparts(obj.files{i});
      end
      udirs=unique(dirs(args.files));

      tref=t(:,:,ref);tref=tref(:);
      tsel=tref>=args.refrange(1) & tref<=args.refrange(2);
      trefs=sort(tref);
      setfig('checktime');clf;
      tiledlayout('flow');
      map=[];
      for j=1:length(udirs);
        nexttile;
        h=[];leg={};
        for i=args.files
          if i==ref
            continue;
          end
          if strcmp(dirs{i},udirs{j})
            t2=t(:,:,i);t2=t2(:);
            if sum(isfinite(t2(tsel)))<4
              continue;
            end
            fit=piecewise(tref(tsel),t2(tsel),args.timetol,2);
            pred=interp1(fit(:,1),fit(:,2),trefs,'linear','extrap');
            fprintf('%-20.20s  [%s]\n',obj.samples{i}, sprintf('(%.2f@%.2f) ',fit'));
            h(end+1)=plot(tref,t2-tref,'o');
            hold on;
            plot(trefs,pred-trefs,'-','Color',get(h(end),'Color'));
            leg{end+1}=obj.samples{i};
          end
        end
        sel=strcmp(dirs,udirs{j});
        tmed=nanmedian(t(:,:,sel),3);
        fit=piecewise(tref(tsel),tmed(tsel),args.timetol,2);
        [~,dirname]=fileparts(udirs{j});
        fprintf('%-20.20s  [%s] over %s\n\n', '', sprintf('(%.2f@%.2f) ',fit'),dirname);
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

      h1=[];h2=[];
      sens=[];
      for i=1:length(obj.files)
        sel=all(obj.contains(:,[i,ref]),2);
        if any(sel)
          sens(i)=nanmedian(reshape(obj.ic(sel,:,i)./obj.ic(sel,:,ref),[],1));
        else
          sens(i)=1;
        end
      end
      ti=['File Sensitivity'];
      setfig(ti);clf;
      bar(sens);
      obj.fsens=sens;
      set(gca,'YScale','log');
      set(gca,'XTick',1:length(obj.samples));
      set(gca,'XTickLabel',obj.samples);
      set(gca,'XTickLabelRotation',90);
      ylabel('Sensitivity');
      title(ti);

      for k=1:length(obj.ADDUCTS)
        % Check sensitivity by target
        for i=1:length(obj.names)
          sel=obj.contains(i,:);
          tsens(i)=nanmedian(squeeze(obj.ic(i,k,sel))./sens(sel)');
        end
        ti=['Target Sensitivity - ',obj.ADDUCTS(k).name];
        setfig(ti);clf;
        bar(tsens);
        h2(k)=gca;
        set(gca,'YScale','log');
        ticks=1:80:length(tsens);
        set(gca,'XTick',ticks)
        set(gca,'XTickLabel',obj.names(ticks));
        set(gca,'XTickLabelRotation',90);
        ylabel('Sensitivity (Ion count for ref file)');
        title(ti);
        ax(:,k)=axis;
        for i=1:size(obj.normic,1)
          obj.normic(i,k,:)=obj.ic(i,k,:)/tsens(i);
        end
        for i=1:size(obj.normic,3)
          obj.normic(:,k,i)=obj.normic(:,k,i)/sens(i);
        end
        obj.tsens(:,k)=tsens;
      end
      linkaxes(h2);

      % Show relative sensitivity of targets
      setfig('Rel Target Sensitivity');clf;
      rel=obj.tsens;
      total=sum(rel,2);
      for i=1:size(rel,2)
        rel(:,i)=rel(:,i)./total;
      end
      bar(rel,'stacked');
      set(gca,'XTick',ticks)
      set(gca,'XTickLabel',obj.names(ticks));
      set(gca,'XTickLabelRotation',90);
      legend({obj.ADDUCTS.name});

      % Cross-scatter of adduct sensitivities
      setfig('Ion Count compare');clf;
      t=tiledlayout('flow');
      t.Title.String='Ion Count Compare';
      h=[];
      for k=2:length(obj.ADDUCTS)
        nexttile;
        h(end+1)=loglog(obj.tsens(:,1),obj.tsens(:,k),'.','MarkerSize',8);
        hold on;
        ax=axis;
        plot(ax(1:2),ax(1:2),':');
        xlabel(sprintf('Ion Count %s',obj.ADDUCTS(1).name));
        ylabel(sprintf('Ion Count %s',obj.ADDUCTS(k).name));
      end
      
      % Overall discrimination
      setfig('Discrimination');clf;
      sel=isfinite(obj.meantime);
      boxplot(reshape(max(obj.normic(sel,:,:),[],2),1,[]),reshape(obj.contains(sel,:),1,[]))
      set(gca,'YScale','log');
      logticks(0,1);
      ylabel('Normalized Ion Count');
      set(gca,'XTickLabel',{'Unexpected','Expected'});
    end
    
    function checknormic(obj)
    % Check pattern of normalized ion counts over assigned compounds
      
      for i=1:length(obj.names)
        n=obj.normic(i,1,obj.contains(i,:));
        if all(isfinite(n))
          fprintf('%3d: %.2f-%.2f %6.2f %s\n', i, min(n), max(n), max(n)/min(n), sprintf('%.2f ',n));
        end
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

    function plotmap(obj,varargin)
      defaults=struct('adduct',1);
      args=processargs(defaults,varargin);

    % Plot map of compounds in mass vs time space
      setfig(sprintf('Compound Map [%s]',obj.ADDUCTS(args.adduct).name));clf;
      etime=[];
      for i=1:length(obj.mass)
        etime(i,1)=nanmin(obj.time(i,args.adduct,:));
        etime(i,2)=nanmax(obj.time(i,args.adduct,:));
        mz(i,1)=nanmin(obj.mz(i,args.adduct,:));
        mz(i,2)=nanmax(obj.mz(i,args.adduct,:));
      end
      plot(obj.mztarget([],args.adduct),obj.meantime,'.b');
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
      title(sprintf('Compound Map [%s] (%d distinguishable/%d identified)',obj.ADDUCTS(args.adduct).name,ngood,sum(isfinite(obj.meantime))));
    end
    
    function ploteics(obj,name,varargin)
    % Plot EIC's for given compound using provided ms cell array aligned with obj.samples
      defaults=struct('mzdata',[],'adduct',1,'falsethresh',0.1,'minic',400,'zoom',true,'mztol',.01,'timerange',[-inf,inf]);
      args=processargs(defaults,varargin);

      if ischar(name)
        ind = obj.find(name);
      else
        ind=name;
      end

      if args.zoom & isfinite(obj.meantime(ind))
        args.timerange=obj.meantime(ind)+[-5,5];
      end
        
      sel=isfinite(squeeze(obj.ic(ind,args.adduct,:)));   % Files that have hits
      contains=obj.contains(ind,:)';

      tp=sel&contains;
      fp=sel&~contains;
      fn=~sel&contains;
      fprintf('Have %d true hits, %d false positives, %d false negatives\n', sum(tp), sum(fp), sum(fn));
      tolist=[find(tp);find(fn);find(fp)];
      ti=sprintf('%s[%s] m/z=%.4f T=%.2f',obj.names{ind},obj.ADDUCTS(args.adduct).name, obj.mztarget(ind,args.adduct), obj.meantime(ind));
      setfig(ti);clf;
      t=tiledlayout('flow');
      t.TileSpacing = 'compact';
      t.Padding = 'compact';
      maxic=max(obj.ic(ind,args.adduct,tolist));

      h=[];
      for ii=1:length(tolist)
        i=tolist(ii);
        nexttile;
        
        filetime=interp1(obj.maps{i}.time(:,1),obj.maps{i}.time(:,2),obj.meantime(ind),'linear','extrap');
        filetimerange=interp1(obj.maps{i}.time(:,1),obj.maps{i}.time(:,2),args.timerange,'linear','extrap');
        filemz=interp1(obj.maps{i}.mz(:,1),obj.maps{i}.mz(:,2),obj.mztarget(ind,args.adduct),'linear','extrap');
        plot(filetime,obj.tsens(ind,args.adduct)*obj.fsens(i),'*','HandleVisibility','off');
        hold on;

        if isempty(args.mzdata) || isempty(args.mzdata{i})
          obj.allfeatures(i).ploteic(filemz,'mztol',args.mztol,'timerange',filetimerange);
        else
          args.mzdata{i}.ploteic(filemz,'newfig',false,'mztol',args.mztol,'timerange',filetimerange);
        end
        h(end+1)=gca;
        xlabel('');ylabel('');
        if ~contains(i) & sel(i)
          title([obj.samples{i},' (False Positive)']);
        elseif contains(i) & ~sel(i)
          title([obj.samples{i},' (False Negative)']);
        else
          title(obj.samples{i});
        end
      end
      linkaxes(h);
      xlabel(t,'Time');
      ylabel(t,'Intensity');
      title(t,ti);
      if args.zoom & isfinite(obj.meantime(ind))
        % Zoom to RT of interest
        ax=axis;
        ax([1,2])=[-2,2]+obj.meantime(ind);
        axis(ax);
      end
    end
    
    function s=featstring(obj,ind,j,k,fi)
    % Check if this fi is assigned to a compound
      [a,b,c]=ind2sub(size(obj.featureindex(:,:,j)),find(obj.featureindex(:,:,j)==fi));
      if isempty(a)
        label='';
      else
        label=obj.names{a};
      end
      s=sprintf(' [%-4d %s %s]',fi,obj.reffeatures(j).features(fi).tostring('mztarget',obj.mztarget(ind,k),'timetarget',obj.time(ind,k),'intensitytarget',obj.fsens(j)*obj.tsens(ind,k),'details',false,'fixedwidth',true),label);
    end
    
    function getinfo(obj,name,varargin)
      defaults=struct('mzdata',[],'adduct',[],'falsethresh',0.1,'minic',400);
      args=processargs(defaults,varargin);

      if ischar(name)
        ind = obj.find(name);
      else
        ind=name;
      end
      if isempty(args.adduct)
        % Check which adduct was used to label, if any
        if ~isempty(obj.astats) && ~isempty(obj.astats(ind).adduct)
          args.adduct=obj.astats(ind).adduct;
        else
          % Show all adducts
          for i=1:length(obj.ADDUCTS)
            obj.getinfo(name,'adduct',i,'mzdata',args.mzdata);
            fprintf('\n');
          end
          return;
        end
      end
      k=args.adduct;
      
      meanic=nanmean(obj.ic(ind,k,obj.contains(ind,:)));
      minic=nanmin(obj.normic(ind,k,obj.contains(ind,:)));
      fprintf('%s[%s] (%d): m/z=%8.4f t=%.2f [%.2f-%.2f] sens=%.0f',obj.names{ind},obj.ADDUCTS(k).name, ind, obj.mztarget(ind,k),obj.meantime(ind),obj.timewindow(ind,:),obj.tsens(ind,k));
      if ~isempty(obj.astats) && ~isempty(obj.astats(ind).run) && obj.astats(ind).adduct==args.adduct
        s=obj.astats(ind);
        fprintf(', FP=%d, FN=%d, hits=%dg,%dL,%dH, miss=%ds,%dw',s.FP,s.FN,s.hitgood,s.hitlow,s.hithigh,s.missstrong,s.missweak);
      end
      fprintf('\n');
      aliases=[];
      label=true;
      for k1=1:length(obj.ADDUCTS)
        aliases=setdiff(find(abs(obj.mass+obj.ADDUCTS(k1).mass-(obj.mass(ind)+obj.ADDUCTS(k).mass))<obj.MZFUZZ*2),ind);
        if length(aliases)>0
          if label
            fprintf('Aliases:\n');
            label=false;
          end
          for ii=1:length(aliases)
            i=aliases(ii);
            imeanic=nanmean(obj.ic(i,obj.contains(i,:)));
            fprintf('\t%-6.6s[%-4.4s]: m/z=%8.4f (d=%4.0f) t=%.2f (d=%.2f) meanic=%.0f\n',...
                    obj.names{i},obj.ADDUCTS(k1).name, obj.mass(i)+obj.ADDUCTS(k1).mass,...
                    (obj.mass(i)+obj.ADDUCTS(k1).mass-(obj.mass(ind)+obj.ADDUCTS(k).mass))*1e4,...
                    obj.meantime(i),obj.meantime(i)-obj.meantime(ind),imeanic);
          end
        end
      end
      if ~isempty(args.mzdata)
        setfig([obj.names{ind},'-',obj.ADDUCTS(args.adduct).name]);
        t=tiledlayout('flow');
        title(t,sprintf('%s m/z=%.4f t=%.2f',obj.names{ind},obj.mztarget(ind,k),obj.meantime(ind)));
      end
      
      for j=1:length(obj.files)
        if ~obj.contains(ind,j)
          continue;
        end
        fprintf('%-15.15s %5.2f',obj.samples{j},obj.fsens(j));

        fi=obj.featureindex(ind,k,j);
        if isfinite(fi) && fi>0
          fprintf('%s',obj.featstring(ind,j,k,fi));
        else
          fprintf(' Expected IC=%6.0f%s',obj.tsens(ind,k)*obj.fsens(j),blanks(40));  % To pad so others align
        end
        others=setdiff(obj.multihits{ind,k,j},fi);
        if length(others)>0
            fprintf('   Others:');
        end
        for p=1:length(others)
          feat=obj.reffeatures(j).features(others(p));
          fprintf('%s',obj.featstring(ind,j,k,others(p)));
        end
        if ~isempty(args.mzdata)
          nexttile;
          obj.plotscan(ind,args.mzdata{j},'adduct',args.adduct);
        end
        fprintf('\n');
      end
      if isfinite(meanic) && isfinite(obj.meantime(ind))
        % False positives
        firstfalse=true;
        for j=1:length(obj.files)
          fi=obj.featureindex(ind,k,j);
          if isfinite(fi) && ~obj.contains(ind,j) && obj.normic(ind,k,j)>=args.falsethresh
            if firstfalse
              fprintf('False positives with normIC >= %.3f:\n',args.falsethresh);
              firstfalse=false;
            end
            fprintf(' %-14.14s %5.2f%s\n',obj.samples{j},obj.fsens(j),obj.featstring(ind,j,k,fi));
            if ~isempty(args.mzdata)
              nexttile;
              obj.plotscan(ind,args.mzdata{j},'adduct',args.adduct);
              ti=get(gca,'Title');
              ti.String=[ti.String,'[Unexp]'];
              set(gca,'Title',ti);
            end
          end
        end
      end
    end

    function plotscan(obj,ind,mzdata,varargin)
      defaults=struct('adduct',1);
      args=processargs(defaults,varargin);

      meanic=nanmean(obj.ic(ind,obj.contains(ind,:)));
      [ic,mz,t]=mzdata.mzscan(obj.mztarget(ind,args.adduct),'mztol',obj.MZFUZZ);
      %plot(t,ic,'o-','MarkerSize',2);
      stem(t,ic,'-','MarkerSize',2);
      hold on;
      set(gca,'YScale','log');
      ax=axis;
      ax(3)=100; ax(4)=max(ax(4),1e5); axis(ax);
      plot(obj.timewindow(ind,1)*[1,1],ax(3:4),':b');
      plot(obj.timewindow(ind,2)*[1,1],ax(3:4),':b');
      if all(isfinite(obj.timewindow(ind,:)))
        set(gca,'XTick',round(obj.timewindow(ind,:)));
      end
      ylabel('Ion Count');
      yyaxis right
      plot(t,mz,'ro','MarkerSize',2);
      hold on;
      ax=axis;
      axis(ax);
      mzrange=obj.mztarget(ind,args.adduct)+obj.MZFUZZ*[-1,1];
      plot(ax(1:2),mzrange(1)*[1,1],':r');
      plot(ax(1:2),mzrange(2)*[1,1],':r');
      ylabel('M/Z');
      if all(isfinite(obj.timewindow(ind,:)))
        ax(1:2)=[obj.timewindow(ind,1)-obj.TIMEFUZZ/2,obj.timewindow(ind,2)+obj.TIMEFUZZ/2];
      end
      ax(3:4)=obj.mztarget(ind,args.adduct)+obj.MZFUZZ*2*[-1,1];
      axis(ax);
      set(gca,'YTick',mzrange);
      title(sprintf('%s - %s',mzdata.name,obj.ADDUCTS(args.adduct).name));
    end
    
    function overlaps(obj)
      [stime,ord]=sort(obj.timewindow(:,1));
      nover=0;
      for i=1:length(stime)
        if isnan(stime(i))
          break;
        end
        ii=ord(i);
        for j=i+1:length(stime)
          if isnan(stime(j))
            break;
          end
          jj=ord(j);
          if obj.timewindow(jj,1)>=obj.timewindow(ii,2)
            break;
          end
          assert(obj.timewindow(ii,2)>obj.timewindow(jj,1));
          for k1=1:length(obj.ADDUCTS)
            m1=(obj.mass(ii)+obj.ADDUCTS(k1).mass);
            for k2=1:length(obj.ADDUCTS)
              m2=(obj.mass(jj)+obj.ADDUCTS(k2).mass);
              mdiff=m1-m2;
              if abs(mdiff) < 2*obj.MZFUZZ
                fprintf('Overlap: %6.6s+%-4.4s (%.4f,[%.2f,%.2f]) and %6.6s+%-4.4s (%.4f,[%.2f,%.2f]) dmz=%.4f\n', obj.names{ii},obj.ADDUCTS(k1).name,m1,obj.timewindow(ii,:),obj.names{jj},obj.ADDUCTS(k2).name,m2,obj.timewindow(jj,:),abs(mdiff));
                nover=nover+1;
              end
            end
          end
        end
      end
      fprintf('Have %d overlaps (%.3f/compound)\n', nover, nover*2/length(obj.names));
    end
    
    function listunexpected(obj,varargin)
    % List peaks that shouldn't be present, in order of descending ion count
      defaults=struct('nlist',30,'adduct',1);
      args=processargs(defaults,varargin);

      fprintf('Comparing for adduct %s:\n', obj.ADDUCTS(args.adduct).name);
      ic=obj.ic;
      ic(obj.contains)=0;
      normic=squeeze(obj.normic(:,args.adduct,:));
      normic(obj.contains)=0;
      [~,ord]=sort(normic(:),'desc','MissingPlacement','last');
      mztarget=obj.mztarget;
      mztarget=mztarget(:,args.adduct);
      for i=1:args.nlist
        [ii,ij]=ind2sub(size(obj.contains),ord(i));
        fprintf('%12.12s:%-7.7s IC=%8.0f(%8.3f)', obj.samples{ij}, obj.names{ii}, ic(ii,args.adduct,ij),normic(ord(i)));
        alias=find(obj.contains(:,ij) & abs(mztarget-mztarget(ii))<=obj.MZFUZZ & abs(obj.meantime-obj.meantime(ii))<=obj.TIMEFUZZ);
        if ~isempty(alias)
          fprintf(' Indistinguishable from ');
          for j=1:length(alias)
            fprintf(' %s',obj.names{alias(j)});
          end
        end
        fprintf('\n');
      end
    end
    
    function conflated(obj,varargin)
    % List possibly conflated compounds
      defaults=struct('adduct',1,'mztol',obj.MZFUZZ,'timetol',obj.TIMEFUZZ);
      args=processargs(defaults,varargin);

      fprintf('mztol=%.4f, timetol=%.2f\n', args.mztol, args.timetol);
      for i=1:length(obj.names)
        if isnan(obj.meantime(i))
          continue;
        end
        for j=i+1:length(obj.names)
          if isnan(obj.meantime(j))
            continue;
          end
          toverlap=abs(obj.meantime(i)-obj.meantime(j)') < args.timetol;
          if toverlap
            mzoverlap=abs(obj.mztarget(i)-obj.mztarget(j)') < args.mztol;
            if any(mzoverlap(:))
              [a2,a1]=ind2sub(size(mzoverlap),find(mzoverlap,1));
              fprintf('%6s[%-4s] %.4f @ %5.2f IC=%6.0f  overlaps %6s[%-4s] %.4f @ %5.2f IC=%6.0f ', ...
                      obj.names{i}, obj.ADDUCTS(a1).name, obj.mztarget(i,a1),obj.meantime(i),obj.tsens(i,a1),...
                      obj.names{j}, obj.ADDUCTS(a2).name, obj.mztarget(j,a1),obj.meantime(j),obj.tsens(j,a2));
              fprintf('dmz=%.4f, dT=%.2f\n',  abs(obj.mztarget(i,a1)-obj.mztarget(j,a2)),abs(obj.meantime(i)-obj.meantime(j)));
            end
          end
        end
      end
      
      nprint=0;
      for k=1:length(obj.samples)
        f=obj.featureindex(:,:,k);  f=f(isfinite(f(:)));f=f(f~=0);
        fs=sort(f);
        dupes=fs(find(diff(fs)==0));
        for i=1:length(dupes)
          [t1,a1]=ind2sub(size(obj.featureindex(:,:,k)),find(obj.featureindex(:,:,k)==dupes(i)));
          fprintf('In %s, feature %d is shared by:\n', obj.samples{k}, dupes(i));
          assert(length(t1)<10);
          for j=1:length(t1)
            fprintf('  %s[%s](%d) mz=%.4f, T=%.2f\n',obj.names{t1(j)}, obj.ADDUCTS(a1(j)).name, t1(j), obj.mztarget(t1(j),a1(j)), obj.meantime(t1(j)));
            nprint=nprint+1;
          end
        end
        if nprint>100
          fprintf('stopping after 100 items\n');
          break;
        end
      end
      
    end
    
    function f=getformula(obj,i)
      f=obj.sdf.getformula(i);
    end

    function export2mzmine(obj,filename,varargin)
    % Export to mzmine2 via CSV
    % Each row has fields:  ID, m/z, retention time, identity, formula
      defaults=struct('adducts',[],'notimes',false,'ind',[]);
      args=processargs(defaults,varargin);
    
      if isempty(args.ind)
        ind=true(size(obj.mass));
      end
      if islogical(args.ind)
        args.ind=find(args.ind);
      end
      if isempty(args.adducts)
        args.adducts=1:length(obj.ADDUCTS);
      end
      fd=fopen(filename,'w');
      fprintf(fd,'ID, m/z, retention time, identity, formula\n');
      for ii=1:length(args.ind)
        i=args.ind(ii);
        for j=1:length(args.adducts)
          t=obj.meantime(i);
          if isnan(t) || args.notimes
            t=0;
          end
          fprintf(fd,'%d,%.4f,%.2f,%s[%s],%s\n', i, obj.mass(i)+obj.ADDUCTS(j).mass,t/60.0,obj.names{i},obj.ADDUCTS(j).name,obj.getformula(i));
        end
      end
      fclose(fd);
    end
    
  end
end
