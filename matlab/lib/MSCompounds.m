% Data structure to hold information about compounds located in mass spec runs
classdef MSCompounds < handle
  properties
    compound; % compound(i) - pk of compound i in compounds database
    names;   % names{i} - Name of compound i (e.g.'51B07')
    mass;    % mass(i) - monoisotpic mass of compound i
    formula; % formula{i} - molecular formula of compound i
    samples;   % samples{j} - Name of sample (mass spec file)
    files;   % files{j} - Mass spec filename j
    maps;    % maps{j} - Piecewise linear maps for converting between references (col 1) and file values (col2) in file j
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
    MZFUZZ;   % global mztol
    TIMEFUZZ; % global timetol
    ADDUCTS;  % adducts to use
  end
  
  properties(Constant)
    dbhost='35.203.151.202';
    dbuser='ngsreadonly';
    dbpassword='';
    database='compounds';
  end
  
  methods
    function obj=MSCompounds(mzfuzz,timefuzz)
      obj.names={};
      obj.formula={};
      obj.multihits={};
      obj.contains=false(0,0);
      obj.samples={};
      obj.files={};
      obj.allfeatures=FeatureList.empty;
      obj.reffeatures=FeatureList.empty;
      obj.astats=struct('run',{},'args',{},'adduct',{},'sel',{},'hitgood',{},'hitlow',{},'hithigh',{},'missstrong',{},'missweak',{},'FP',{},'FN',{},'fpsel',{});
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

    function q=copypart(obj,varargin)
    % Create a new object containing only the selected parts
      defaults=struct('fsel',1:length(obj.files),'csel',1:length(obj.names),'asel',1:length(obj.ADDUCTS));
      args=processargs(defaults,varargin);

      if islogical(args.fsel)
        assert(length(args.fsel)==length(obj.files));
        args.fsel=find(args.fsel);
      else
        assert(all(args.fsel>=1 & args.fsel<=length(obj.files)));
      end
      if islogical(args.csel)
        assert(length(args.csel)==length(obj.names));
        args.csel=find(args.csel);
      else
        assert(all(args.csel>=1 & args.csel<=length(obj.names)));
      end
      if islogical(args.asel)
        assert(length(args.asel)==length(obj.ADDUCTS));
        args.asel=find(args.asel);
      else
        assert(all(args.asel>=1 & args.asel<=length(obj.ADDUCTS)));
      end
      q=MSCompounds(obj.MZFUZZ, obj.TIMEFUZZ);
      q.compound=obj.compound(args.csel);
      q.names=obj.names(args.csel);
      q.mass=obj.mass(args.csel);
      q.formula=obj.formula(args.csel);
      q.samples=obj.samples(args.fsel);
      q.files=obj.files(args.fsel);
      q.maps=obj.maps(args.fsel);
      q.moles=obj.moles(args.fsel);
      q.group=obj.group(args.fsel);
      q.contains=obj.contains(args.csel,args.fsel);
      q.mz=obj.mz(args.csel,args.asel,args.fsel);
      q.time=obj.time(args.csel,args.asel,args.fsel);
      q.meantime=obj.meantime(args.csel);
      q.timewindow=obj.timewindow(args.csel,:);
      q.ic=obj.ic(args.csel,args.asel,args.fsel);
      q.normic=obj.normic(args.csel,args.asel,args.fsel);
      if ~isempty(obj.multihits)
        q.multihits=obj.multihits(args.csel,args.asel,args.fsel);
      end
      if ~isempty(obj.allfeatures)
        q.allfeatures=obj.allfeatures(args.fsel);
      end
      if ~isempty(obj.reffeatures)
        q.reffeatures=obj.reffeatures(args.fsel);
      end
      if ~isempty(obj.featureindex)
        q.featureindex=obj.featureindex(args.csel,args.asel,args.fsel);
      end
      q.tsens=obj.tsens(args.csel,args.asel);
      q.fsens=obj.fsens(args.fsel);
      if ~isempty(obj.astats)
        q.astats=obj.astats(args.csel);
      end
      q.ADDUCTS=obj.ADDUCTS(args.asel);
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
    
    function fl=checkComposition(obj,ms,varargin)  % TODO - test
    % Check composition in mass spec file using current set of compounds
      defaults=struct('debug',false,'timetol',obj.TIMEFUZZ,'mztol',obj.MZFUZZ,'noise',500);
      args=processargs(defaults,varargin);
      sel=find(isfinite(obj.meantime));
      adduct=[obj.astats(sel).adduct];
      mztarget=arrayfun(@(z) obj.mztarget(sel(z),adduct(z)),1:length(sel));
      fl=ms.targetedFeatureDetect(mztarget,'rt',obj.meantime(sel),'mztol',args.mztol,'timetol',args.timetol,'names',obj.names(sel),'noise',args.noise,'debug',args.debug);
      for i=1:length(fl.features)
        f=fl.features(i);
        if length(f.labels)==1
          ind=find(strcmp(f.labels{1},obj.names(sel)));
          assert(length(ind)==1);
          normic=f.intensity/obj.tsens(sel(ind),adduct(ind));
        else
          normic=nan;
        end
        f.extra=struct('normic',normic);
      end
    end
    
    function plotComposition(obj,ms,varargin)
      defaults=struct('debug',false,'timetol',obj.TIMEFUZZ,'mztol',obj.MZFUZZ,'noise',500);  
      args=processargs(defaults,varargin);

      fl=obj.checkComposition(ms,'debug',args.debug,'timetol',args.timetol,'mztol',args.mztol,'noise',args.noise);  
      data=[];
      plate=arrayfun(@(z) str2num(z.name(1:end-3)),fl.features);
      uplates=unique(plate);
      colname={};
      for i=1:length(fl.features)
        f=fl.features(i);
        nm=f.name;
        pnum=find(plate(i)==uplates);
        row=nm(end-2)-'A'+1;
        col=str2num(nm(end-1:end));
        cnum=(row-1)+(col-2)*8+1;
        data(pnum,cnum)=f.extra.normic;
        colname{cnum}=sprintf('%c%02d',row-1+'A',col);
      end
      data(end+1,:)=nan;
      data(:,end+1)=nan;
      ti=['Composition ',fl.name];
      setfig(ti);clf;
      pcolor(log10(data));
      colorbar;
      set(gca,'XTick',(1:8:length(colname))+0.5);
      set(gca,'XTickLabel',colname(1:8:end));
      set(gca,'XTickLabelRotation',90);
      set(gca,'YTick',(1:4:length(uplates))+0.5);
      set(gca,'YTickLabel',uplates(1:4:end));
      title(ti);
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
    
    function [matic,id,refid]=plotCompositionOld(obj,ms,varargin)   % TODO - fix
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
      
    function addCompound(obj,compound,name,mass,formula)
    % Add a compound
      if any(strcmp(obj.names,name)) || (isfinite(compound) && ismember(compound,obj.compound))
        error('%s (%d) already added', name, compound);
      end

      obj.names{end+1}=name;
      nindex=length(obj.names);
      obj.compound(nindex)=compound;
      obj.mass(nindex)=mass;
      obj.formula{nindex}=formula;
      obj.tsens(nindex,1:length(obj.ADDUCTS))=nan;
      obj.astats(nindex)=struct('run',{[]},'args',{[]},'adduct',{[]},'sel',{[]},'hitgood',{[]},'hitlow',{[]},'hithigh',{[]},'missstrong',{[]},'missweak',{[]},'FP',{[]},'FN',{[]},'fpsel',{[]});
      if length(obj.files)>0
        obj.mz(nindex,:,:)=nan;
        obj.time(nindex,:,:)=nan;
        obj.timewindow(nindex,:)=nan;
        obj.meantime(nindex)=nan;
        obj.ic(nindex,:,:)=nan;
        obj.normic(nindex,:,:)=nan;
        obj.contains(nindex,:)=false;
      else
        obj.mz=nan(nindex,length(obj.ADDUCTS), 0);
        obj.time=nan(nindex,length(obj.ADDUCTS), 0);
        obj.timewindow=nan(nindex,0);
        obj.meantime=nan(nindex,1);
        obj.ic=nan(nindex,length(obj.ADDUCTS), 0);
        obj.normic=nan(nindex,length(obj.ADDUCTS), 0);
        obj.contains=false(nindex,0);
      end
    end
    
    function addCompoundsFromSDF(obj,sdf)
    % Add all the compounds in the given SDF file using a name formed from the PLATE and WELL
      assert(isempty(obj.names));   % Otherwise, can't set obj.sdf to correspond
      for i=1:length(sdf.sdf)
        s=sdf.sdf(i);
        name=sprintf('%d%s',str2num(s.BATCH_PLATE(5:end)),s.BATCH_WELL);
        obj.addCompound(nan, name,s.MostAbundantMass,s.getformula());  % May be different from monoisotopic mass
        index=strcmp(obj.names,name);
      end
    end

    function dbopen(data)
    % Open a database connection using given user, close all other connections
      if mysql('status')==0
        res=mysql('select USER();');
        sp=split(res{1},'@');
        if strcmp(sp{1},data.dbuser)
          mysql('use',data.database);
          return;
        end
      end
      mysql('closeall');
      mysql('open',data.dbhost,data.dbuser,data.dbpassword);
      mysql('use',data.database);
    end

    function addCompoundsFromDB(obj,ids)
    % Add the compounds in the given list of db pks (or all if none set)
      obj.dbopen();
      cmd='SELECT compound,name,monoisotopicMass,formula FROM compounds.compounds';
      if nargin>=2
        idlist=sprintf('%d,',ids);
        idlist=idlist(1:end-1);  % Remove trailing comm
        cmd=sprintf('%s WHERE compound IN (%s)',cmd, idlist);
      end
      [compound,name,mass,formula]=mysql(cmd);
      for i=1:length(compound)
        obj.addCompound(compound(i),name{i},mass(i),formula{i});
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
        allid=obj.checkComposition(ms,'map',map,'timetol',500/60/iter,'mztol',obj.MZFUZZ*2/iter);
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
      defaults=struct('debug',0,'timetol',obj.TIMEFUZZ,'minhits',3,'mztol',obj.MZFUZZ,'plot','','minic',1000,'trace',[],'maxFN',0,'maxFP',0,'clear',true,'detectionThreshold',2000,'normicrange',[0.4,2.5],'falsethresh',0.1,'usefiles',[],'allowambig',true);
      args=processargs(defaults,varargin);

      if args.clear || isempty(obj.astats)
        arun=1;
      else
        arun=max([obj.astats.run])+1;
      end
      
      fprintf('assignTimes (run %d):\n',arun);
      if args.clear
        obj.ic(:)=nan;
        obj.normic(:)=nan;
        obj.mz(:)=nan;
        obj.time(:)=nan;
        obj.meantime(:)=nan;
        obj.timewindow(:,:)=nan;
        obj.featureindex=nan(length(obj.names),length(obj.ADDUCTS),length(obj.samples));
        obj.astats(:)=struct('run',{[]},'args',{[]},'adduct',{[]},'sel',{[]},'hitgood',{[]},'hitlow',{[]},'hithigh',{[]},'missstrong',{[]},'missweak',{[]},'FP',{[]},'FN',{[]},'fpsel',{[]})
      end
      used=cell(size(obj.files));
      for j=1:length(obj.files)
        used{j}=obj.featureindex(:,:,j);
        used{j}(~obj.contains(:,j),:)=nan;   % Don't mark unexpected hits as used
        used{j}=used{j}(isfinite(used{j}(:)));
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
              %used=obj.featureindex(:,:,j);
              %used=used(isfinite(used(:)));
              featinds=mhits(~ismember(mhits,used{j}));   % This is faster than using setdiff()
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
                timewindow=mean(esort([m,n]))+args.timetol*[-1,1];
                sel=find(etimes>=timewindow(1) & etimes<=timewindow(2));
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
                expectedIC=tsens*obj.fsens;
                normic=zeros(size(obj.fsens));
                normic(srcfile(sel))=intensity(sel)./expectedIC(srcfile(sel));
                assert(all(isfinite(normic)));
                % Categorize each file as
                %  trueneg - unexpected and not present
                %  falsepos - unexpected but present and > falsethresh
                %  falseweak - unexpected but present and < falsethresh
                %  hitgood - expected and present and in correct normic range
                %  hitlow - expected and present but low normic
                %  hithigh - expected and present but high normic
                %  missweak - expected (weakly) and not present
                %  missstrong - expected (strongly) and not present
                present=ismember(1:length(obj.samples),srcfile(sel));
                trueneg=~present&~expected;
                lowic=normic<args.normicrange(1);
                highic=normic>args.normicrange(2);
                falsepos=present&~expected&normic>=args.falsethresh;
                falseweak=present&~expected&normic<args.falsethresh;
                hit = present & expected;
                hitlow=hit & lowic; 
                hithigh=hit &highic;
                hitgood=hit & ~hitlow & ~hithigh; 
                miss=~present & expected;
                missstrong=miss&~(expectedIC<args.detectionThreshold);   % Include nan expectedIC here
                missweak=miss&expectedIC<args.detectionThreshold;

                assert(all(trueneg|falsepos|falseweak|hitgood|hitlow|hithigh|missweak|missstrong));
                assert(sum(trueneg)+sum(falsepos)+sum(falseweak)+sum(hitgood)+sum(hitlow)+sum(hithigh)+sum(missweak)+sum(missstrong)==length(obj.samples));
                nFP=sum(falsepos);
                nFN=sum(missstrong|hitlow|hithigh);
                if ismember(i,args.trace)
                  selstr='';
                  for ij=1:length(sel)
                    selstr=sprintf('%s %.2f',selstr,etimes(sel(ij)));
                    if any(sel(ij)==selexpected)
                      selstr(end+1)='*';
                    end
                  end
                  if length(selstr)>0
                    selstr=selstr(2:end);
                  end
                  fprintf('m=%d,n=%d, keep=%s, hits=(%d good,%d low, %d high), miss=(%d strong,%d weak), FP=%d,FN=%d, width=%.2f,rng=[%.2f,%.2f]\n',m,n,selstr,sum(hitgood),sum(hitlow),sum(hithigh),sum(missstrong),sum(missweak),nFP,nFN,esort(n)-esort(m),esort([m,n]));
                end
                if sum(hitgood)<args.minhits || nFP>args.maxFP || nFN > args.maxFN
                  continue;
                end

                if length(selexpected)>length(bestset) || (length(selexpected)==length(bestset) && (nFN<astats.FN || (nFN==astats.FN && (nFP<astats.FP || (args.allowambig && nFP==astats.FP && tsens>bestic)))))
                  bestset=selexpected;
                  bestwindow=timewindow;
                  bestic=tsens;
                  astats=struct('run',arun,'args',args,'adduct',k,'sel',srcfile(selexpected),'hitgood',sum(hitgood),'hitlow',sum(hitlow),'hithigh',sum(hithigh),'missstrong',sum(missstrong),'missweak',sum(missweak),'FP',nFP,'FN',nFN,'fpsel',srcfile(setdiff(sel,selexpected)));
                  nbest=1;
                  if ismember(i,args.trace)
                    fprintf('best set: FN=%d,FP=%d,tsens=%.0f:  %s\n', nFN, nFP, tsens, sprintf('%d ',selexpected));
                  end
                elseif length(bestset)==length(selexpected) && any(~ismember(selexpected,bestset))
                  % Different set
                  if ismember(i,args.trace)
                    fprintf('disjoint set: FN=%d,FP=%d,tsens=%.0f:  %s\n', nFN, nFP, tsens, sprintf('%d ',selexpected));
                  end
                  nbest=nbest+1;
                else
                  if ismember(i,args.trace)
                    fprintf('ignored set: %s\n', sprintf('%d ',selexpected));
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
          elseif nbest>1 && ~args.allowambig
            if args.debug  || ismember(i,args.trace)
              fprintf('Have %d equivalent elution times for %s\n', nbest, obj.names{i});
            end   
            continue;
          end

          nFP=length(unique(srcfile(~cont & etimes>=bestwindow(1) & etimes<=bestwindow(2))));
          ntrue=length(unique(srcfile(cont & etimes>=bestwindow(1) & etimes<=bestwindow(2))));
          nFN=sum(obj.contains(i,:))-ntrue;
          if args.debug || strcmp(args.plot,obj.names{i})  || ismember(i,args.trace)
            fprintf('Best has %d true, %d FN and %d FP over [%.2f,%.2f]\n',ntrue,nFN,nFP,bestwindow);
          end

          % Construct consensus view
          obj.meantime(i)=mean(bestwindow);
          obj.timewindow(i,1:2)=obj.meantime(i)+args.timetol*[-1,1];
          obj.astats(i)=astats;
          for kk=1:length(obj.ADDUCTS)
            for j=1:length(obj.files)
              fi=obj.multihits{i,kk,j};
              time=[obj.reffeatures(j).features(fi).time];
              sel=find(time>=obj.timewindow(i,1) & time <= obj.timewindow(i,2));
              fi=fi(sel);
              if ~isempty(fi)
                if length(fi)>1
                  % More than one feature in time window; use one with highest peak height
                  if args.debug || ismember(i,args.trace)
                    fprintf('Have %d ambiguous features for %s[%s] in %s\n', length(fi), obj.names{i},obj.ADDUCTS(k).name,obj.samples{j});
                  end
                  %[~,mind]=min(abs([obj.reffeatures(j).features(fi).time]-obj.meantime(i)));
                  [~,mind]=max([obj.reffeatures(j).features(fi).intensity]);
                  fi=fi(mind);
                end
                obj.featureindex(i,kk,j)=fi;
                if obj.contains(i,j)
                  % Only set used if this feature should be present
                  used{j}(end+1)=fi;
                end
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
    
    function assignSummary(obj)
    % Display information about assignments
      for run=1:max([obj.astats.run])
        fprintf('Run %d',run);
        for adduct=1:length(obj.ADDUCTS)
          n(run,adduct)=sum([obj.astats.run]==run & [obj.astats.adduct]==adduct);
          fprintf(' [%s]=%4d',obj.ADDUCTS(adduct).name, n(run,adduct));
        end
        fprintf(' Total=%4d\n',sum(n(run,:)));
      end
      fprintf('All: ');
      for adduct=1:length(obj.ADDUCTS)
        fprintf(' [%s]=%4d',obj.ADDUCTS(adduct).name, sum(n(:,adduct)));
      end
      fprintf(' Total=%4d\n',sum(n(:)));
    end
    
    function remapfeatures(obj)
    % Create remapped version of file features into reference space
      if length(obj.allfeatures)==length(obj.reffeatures)
        fprintf('Feature lists already remapped\n');
        return;
      end
      fprintf('Remapping %d feature lists...',length(obj.allfeatures));
      for i=1:length(obj.allfeatures)
        if length(obj.reffeatures)<i || isempty(obj.reffeatures(i)) || length(obj.reffeatures(i).features)~=length(obj.allfeatures(i).features)
          fprintf('%d...',i);
          obj.reffeatures(i)=obj.allfeatures(i).maptoref(obj.maps{i});
        end
      end
      fprintf('done\n');
    end
    
    function addMS(obj,ms,mixture,varargin)
    % Add the unique peaks for compounds from a particular M/S run
    % Use prior analyses to figure out the expected elution time for each compound
    % or scan all elution times if the no prior data (keep only if a unique peak is determined)
      defaults=struct('group','','map',struct('mz',[0 0; 1 1 ],'time',[0 0; 1 1 ]),'sample',[]);
      % mzmap(i,2) - piecewise linear map for M/Z; mzmap(:,1) is true M/Z, mzmap(:,2) is for values in ms file
      % timemap(i,2) - piecewise linear map for elution times; timemap(:,1) is "standard" elution time
      args=processargs(defaults,varargin);
      
      if ~isstruct(mixture)
        error('addMS now takes a single mixture struct as the mixture arg');
      end
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
      if isempty(mixture.contents)
        error('Mixture %s (%d) has no contents',mixture.name,mixture.mixture);
      end
      obj.contains(:,findex)=ismember(obj.compound,mixture.contents);

      if ~isempty(ms.featurelists)
        fl=ms.featurelists(end);
        obj.allfeatures(findex)=fl;
      else
        fprintf('Warning: not setting MSCompounds.featurelists\n');
      end
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
        mzlist=[obj.reffeatures(findex).features.mz];
        iclist=[obj.reffeatures(findex).features.intensity];
        % Attempt to locate each one uniquely
        maxic=[];   % Maximum IC per target
        fprintf('%s...',obj.samples{findex});
        for i=1:length(obj.mass)
          if mod(i,1000)==1
            fprintf('%d...',i);
          end
          for k=1:length(obj.ADDUCTS)
            % Use features
            findices=find(abs(mzlist-obj.mztarget(i,k))<args.mztol);
            obj.multihits{i,k,findex}=findices;
            maxic(i,k)=max([0,iclist(findices)]);
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
    
    function x=report(obj,varargin) 
    % Build table on data by compound
    % Each row is a single compound
    % Data by group
      defaults=struct('sort',false,'bygroup',false);
      args=processargs(defaults,varargin);

      if any(isfinite(obj.meantime) & ~any(isfinite(obj.tsens),2))
        error('Run checksensitivity before report');
      end
      
      if args.bygroup
        ugroups=unique(obj.group,'sorted');
      else
        ugroups={'all'};
      end
      if args.sort
        [~,ord]=sort(obj.names);
      else
        ord=1:length(obj.names);
      end
      x=[];row=1;
      for ii=1:length(obj.names)
        i=ord(ii);
        for j=1:length(ugroups)
         if args.bygroup
           files=find(strcmp(obj.group,ugroups{j})& obj.contains(i,:));
         else
           files=find(obj.contains(i,:));
         end
         if ~isempty(files)
           if ~isempty(obj.astats) && i<=length(obj.astats) && ~isempty(obj.astats(i).adduct)
             adduct=obj.astats(i).adduct;  % Same adduct as used to assign it
           else
             [~,adduct]=max(obj.tsens(i,:));   % Use highest adduct
           end
           x(row).name=obj.names{i};
           x(row).mass=round(obj.mass(i),5);
           x(row).adduct=obj.ADDUCTS(adduct).name;
           x(row).mz=round(obj.mztarget(i,adduct),5);
           x(row).time=round(obj.meantime(i),2);
           x(row).sensitivity=round(obj.tsens(i,adduct));
           x(row).FP=obj.astats(i).FP;
           x(row).FN=obj.astats(i).FN;
           for m=1:length(files)
             x(row).file{m}=obj.samples{files(m)};
             x(row).mzoffset(m)=round((obj.mz(i,adduct,files(m))-obj.mztarget(i,adduct))*1e5);
             x(row).rt(m)=round(obj.time(i,adduct,files(m)),2);
             x(row).ic(m)=round(obj.ic(i,adduct,files(m)));
           end
           row=row+1;
         end
        end
      end
      x=struct2table(x);
    end
    
    function platesummary(obj)
    % Summarize stats by original CDIV plate
      plate=cellfun(@(z) z(1:end-3),obj.names,'Unif',false);
      uplate=unique(plate);
      for i=1:length(uplate)
        sel=strcmp(plate,uplate{i});
        nfound(i)=sum(isfinite(obj.meantime(sel)));
        sens(i)=nanmedian(nanmax(obj.tsens(sel,:),[],2));
        fprintf('%s: %d found, %.0f sens\n', uplate{i}, nfound(i), sens(i));
      end
      setfig('platesummary');clf;
      tiledlayout('flow');
      nexttile;
      bar(nfound);
      set(gca,'XTick',1:length(uplate));
      set(gca,'XTickLabel',uplate);
      set(gca,'XTickLabelRotation',90);
      ylabel('Number of compounds found (/80)');
      nexttile
      bar(sens);
      set(gca,'XTick',1:length(uplate));
      set(gca,'XTickLabel',uplate);
      set(gca,'XTickLabelRotation',90);
      ylabel('Median sensitivity');
      set(gca,'YScale','log');
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
      defaults=struct('files',1:length(obj.files),'icplot',false);
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
            x=obj.mztarget;y=obj.mz(:,:,i);ic=obj.ic(:,:,i);
            x(~obj.contains(:,i),:)=nan;   % Only ones that it is supposed to contain
            npts=sum(isfinite(y(:)) & isfinite(x(:)));
            if npts<10
              fprintf('Only %d data points for %s ... skipping\n', npts, obj.files{i});
              continue;
            end
            if args.icplot
              h(end+1)=semilogx(ic(:),1e4*(y(:)-x(:)),'.b','MarkerSize',1);
              hold on;
            else
              err=nanmedian(y(:)-x(:));
              fit=robustfit(x(:),y(:));
              fprintf('%-20.20s  %5.1f [%5.1f, %5.1f] N=%d\n',obj.samples{i}, err*1e4,1e4*(fit(1)+(fit(2)-1)*rng),npts);
              h(end+1)=plot(x(:),1e4*(y(:)-x(:)),'o');
              hold on;
              plot(rng,1e4*(fit(1)+(fit(2)-1)*rng),'-','Color',get(h(end),'Color'));
            end
            leg{end+1}=obj.samples{i};
            allx=[allx;x(:)];
            ally=[ally;y(:)];
          end
        end
        if args.icplot
          xlabel('Ion Count');
        else
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
        end
          
        ylabel('File-True M/Z *1e4');
        title(udirs{j});
      end
    end
    
    function [allmap,filemap]=checktime(obj,varargin)
    % Check all multiple hits of ref against each other
      defaults=struct('debug',false,'timetol',obj.TIMEFUZZ,'minic',1e4,'ref',[],'refrange',[5,30],'files',1:length(obj.files),'assign','','current',true);
      args=processargs(defaults,varargin);

      if ~ismember(args.assign,{'','dir','file'})
        error('assign value of ''%s'' is not one of '''',''dir'', or ''file''',args.assign);
      end
      if ~isempty(args.assign) && args.current && any(cellfun(@(z) any(z.time(:,1)~=z.time(:,2)) ,obj.maps))
        error('maps already non-identity and requested to reassign -- need to set current to false');
      end
      if args.current && isempty(obj.reffeatures)
        error('Need to run remapfeatures');
      end
      args.files=union(args.files,args.ref);
      t = nan(length(obj.mass),length(obj.samples));
      ic = t;
      for i=1:length(obj.names)
        if isfinite(obj.meantime(i))
          k=obj.astats(i).adduct;
          for jj=1:length(args.files)
            j=args.files(jj);
            if obj.contains(i,j) && obj.ic(i,k,j)>=args.minic
              fi=obj.featureindex(i,k,j);
              if args.current
                t(i,j)=obj.reffeatures(j).features(fi).time;
              else
                t(i,j)=obj.allfeatures(j).features(fi).time;
              end
            end
          end
        end
      end
      % Setup linear matrix equation to solve for linear time mapping
      nf=length(obj.samples);
      A=zeros(0,2*nf);
      for i=1:length(obj.names)
        for j1=1:size(t,2)
          if isnan(t(i,j1))
            continue;
          end
          for j2=j1+1:size(t,2)
            if isnan(t(i,j2))
              continue;
            end
            A(end+1,[j1,j1+nf,j2,j2+nf])=[t(i,j1),1,-t(i,j2),-1];
          end
        end
      end
      b=zeros(size(A,1),1);
      % Add equations for refs
      if isempty(args.ref)
        % No ref, add equation such that mean slope,intercept is 1,0
        wt=1e6;   % High enough to be sure this dominates
        A(end+1,1:nf)=wt; b(end+1)=nf*wt;
        A(end+1,nf+1:end)=wt; b(end+1)=0;
      else
        for jj=1:length(args.ref)
          j=args.ref(jj);
          A(end+1,j)=1e6; b(end+1)=1e6;
          A(end+1,j+nf)=1; b(end+1)=0;
        end
      end
      asel=find(any(A));
      A=A(:,asel);   % Only columns that we used
      nf=length(asel)/2;
      fprintf('Solving %d equations in %d unknowns with %d refs\n', size(A),length(args.ref));
      x=A\b;

      % m,k are slope intercept;  nan for files not in args.files
      m=nan(length(obj.samples),1);
      k=nan(length(obj.samples),1);
      m(asel(1:nf))=x(1:nf);
      k(asel(1:nf))=x(nf+1:end);
      % Compute tadj (times corrected by linear map)
      tadj=nan(size(t));
      for i=1:length(m)
        tadj(:,i)=t(:,i)*m(i)+k(i);
      end
      
      % Group by source folder (Mass Spec run)
      dirs={};
      for ii=1:length(args.files)
        i=args.files(ii);
        dirs{i}=fileparts(obj.files{i});
      end
      udirs=unique(dirs(args.files));
      
      % Overall map for files
      nd=length(udirs);
      AT=zeros(size(A,1),2*nd);
      for i=1:nd
        sel=find(strcmp(dirs,udirs{i}));
        AT(:,i)=nansum(A(:,sel),2);
        AT(:,i+nd)=nansum(A(:,sel+nf),2);
      end
      xT=AT\b;
      mT=xT(1:nd);
      kT=xT(nd+1:end);
      
      setfig('checktime');clf;
      layout=tiledlayout('flow');
      if args.current
        atype='showing additional mapping needed';
      else
        atype='cumulative mapping';
      end
        
      title(layout,sprintf('Time Alignment N=%d, %s', size(A,1), atype));
      allmap=[];filemap=[];
      
      allgca=[];
      for j=1:length(udirs);
        nexttile;
        allgca(end+1)=gca;
        h=[];leg={};
        mapT=[args.refrange;(args.refrange-kT(j))/mT(j)]';
        filemap=[filemap,struct('dir',udirs{j},'map',mapT)];

        for ii=1:length(args.files)
          i=args.files(ii);
          if strcmp(dirs{i},udirs{j})
            t1=nanmean(tadj(:,[1:i-1,i+1:end]),2);   % Using all others with slope adjustment
            t2=t(:,i);   % This one, unadjusted so we can see the original errors
            pred=tadj(:,i);   % Adjusted
            if sum(isfinite(t2))<4
              continue;
            end
            h(end+1)=plot(t1,t2-t1,'.');   % Error from average of all of the other adjusted ones
            hold on;
            map=[args.refrange;(args.refrange-k(i))/m(i)]';   % map(1,:) is ref, map(2,:) is file values
            fprintf('%-20.20s  [%s] %.3ft+%6.3f\n',obj.samples{i}, sprintf('%.2f->%.2f ',map'),m(i),k(i));
            plot(map(:,1),map(:,2)-map(:,1),'-','Color',get(h(end),'Color'));
            leg{end+1}=obj.samples{i};
            allmap=[allmap,struct('file',obj.files{i},'map',map)];
            if strcmp(args.assign,'file')
              obj.maps{i}.time=map;
            elseif strcmp(args.assign,'dir')
              obj.maps{i}.time=mapT;
            end
          end
        end

        % Plot directory-wide map
        [~,dirname]=fileparts(udirs{j});
        fprintf('%-20.20s  [%s] %.3ft+%6.3f over %s\n\n', '', sprintf('%.2f->%.2f ',mapT'),mT(j), kT(j), dirname);
        h(end+1)=plot(mapT(:,1),mapT(:,2)-mapT(:,1),'k','linewidth',2);
        leg{end+1}='All';

        % Complete plot setup
        plot(args.refrange,obj.TIMEFUZZ*[1,1],'k:','linewidth',1,'HandleVisibility','off');
        plot(args.refrange,-obj.TIMEFUZZ*[1,1],'k:','linewidth',1,'HandleVisibility','off');
        ax=axis;
        ax(3)=nanmin(t(:)-tadj(:)-obj.TIMEFUZZ*2);
        ax(4)=nanmax(t(:)-tadj(:)+obj.TIMEFUZZ*2);
        axis(ax);
        legend(h,leg,'location','best');
        xlabel('Reference time');
        ylabel('File-Ref time');
        title(udirs{j});
      end
      linkaxes(allgca);
      if ~strcmp(args.assign,'')
        fprintf('Time maps were reset; clearing reffeatures -- need to rerun findfeatures');
        obj.reffeatures=FeatureList.empty;
      end
    end
    
    function checksensitivity(obj,varargin)
    % Check sensitivity by file relative to ref file
      defaults=struct('ref',[],'plot',false);
      args=processargs(defaults,varargin);

      obj.normic=obj.ic;

      sens=ones(size(obj.files));
      for pass=1:10
        ratio=nan(length(obj.files));
        for i=1:length(obj.files)
          for j=1:length(obj.files)
            sel=all(obj.contains(:,[i,j]),2);
            if sum(sel)>=3
              ratio(i,j)=nanmedian(reshape(obj.ic(sel,:,i)./obj.ic(sel,:,j),[],1))*sens(j)/sens(i);
            end
          end
        end
        sens=nanmedian(ratio,2)'.*sens;
      end
      if any(isnan(sens))
        fprintf('Setting sensitivity of %d files without enough data to 1.0\n', sum(isnan(sens)));
        sens(isnan(sens))=1;
      end
      
      if ~isempty(args.ref)
        sens=sens/sens(args.ref);
      end
      obj.fsens=sens;

      if args.plot
        ti=['File Sensitivity'];
        setfig(ti);clf;
        tiledlayout('flow');
        nexttile;
        bar(sens);
        set(gca,'YScale','log');
        set(gca,'XTick',1:length(obj.samples));
        set(gca,'XTickLabel',obj.samples);
        set(gca,'XTickLabelRotation',90);
        ylabel('Sensitivity');
        title(ti);
        nexttile;
        boxplot(sens,obj.moles);
        xlabel('Moles injected');
        ylabel('Sensitivity');
      end
      
      h2=[];
      for k=1:length(obj.ADDUCTS)
        % Check sensitivity by target
        for i=1:length(obj.names)
          sel=obj.contains(i,:);
          tsens(i)=nanmedian(squeeze(obj.ic(i,k,sel))./sens(sel)');
        end
        if args.plot
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
        end
        
        for i=1:size(obj.normic,1)
          obj.normic(i,k,:)=obj.ic(i,k,:)/tsens(i);
        end
        for i=1:size(obj.normic,3)
          obj.normic(:,k,i)=obj.normic(:,k,i)/sens(i);
        end
        obj.tsens(:,k)=tsens;
      end
      if args.plot
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
    % Plot map of compounds in mass vs time space
      defaults=struct();
      args=processargs(defaults,varargin);

      setfig('Compound Map');clf;
      etime=obj.timewindow;
      mz=nan(length(obj.mass),2);
      for i=1:length(obj.mass)
        if isfinite(obj.meantime(i))
          mz(i,:)=obj.mztarget(i,obj.astats(i).adduct)+obj.MZFUZZ*[-1,1];
        end
      end
      plot(mean(mz,2),obj.meantime,'.b');
      hold on;
      mzmissing=obj.mass(isnan(obj.meantime));
      plot(mzmissing,-2+2*rand(size(mzmissing)),'.r','MarkerSize',1);
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
      title(sprintf('Compound Map (%d distinguishable/%d identified)',ngood,sum(isfinite(obj.meantime))));
    end
    
    function ploticreplicates(obj,varargin)
    % Compare ion count from each file with average from others
      defaults=struct('fsel',1:length(obj.samples),'outlierthresh',20);
      args=processargs(defaults,varargin);

      data=nan(length(obj.mass),length(obj.samples));
      for i=1:length(obj.mass)
        if isfinite(obj.meantime(i))
          adduct=obj.astats(i).adduct;
          data(i,:)=squeeze(obj.ic(i,adduct,:))./obj.fsens';
        end
      end
      data(~obj.contains)=nan;   % Ony consider ones where the target is supposed to be present
      setfig('IC Replicates');clf;
      fprintf('Conditions with >%.0fx error:\n', args.outlierthresh);
      all=[];
      for jj=1:length(args.fsel)
        j=args.fsel(jj);
        icexp=nanmean(data(:,[1:j-1,j+1:end]),2)*obj.fsens(j);
        loglog(icexp,data(:,j),'.b','MarkerSize',1);
        all=[all;icexp,data(:,j)];
        hold on;
        bad=abs(log10(data(:,j)./icexp))>log10(args.outlierthresh);
        if sum(bad)>0
          loglog(icexp(bad),data(bad,j),'.r');
          fprintf('Sample %2d %s: %s\n', j, obj.samples{j}, sprintf('%4d ',find(bad)));
        end
      end
      avgerr=exp(nanmean(abs(log(all(:,1)./all(:,2)))));
      fprintf('Error sigma = %.2fx\n', avgerr);
      hold on;
      ax=axis;
      plot(ax(1:2),ax(1:2),'r:');
      xlabel('Expected ion count based on all other files');
      ylabel('Normalized ion count from single file');
      title('Ion Count Replicates');
    end
      
    function ploteics(obj,name,varargin)
    % Plot EIC's for given compound using provided ms cell array aligned with obj.samples
      defaults=struct('mzdata',[],'adduct',[],'falsethresh',0.1,'minic',400,'zoom',true,'mztol',obj.MZFUZZ,'timerange',[-inf,inf]);
      args=processargs(defaults,varargin);

      if ischar(name)
        ind = obj.find(name);
      else
        ind=name;
      end

      if isempty(args.adduct)
        args.adduct=obj.astats(ind).adduct;
        if isempty(args.adduct)
          args.adduct=1;
        end
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
      s=sprintf(' [%-5d %s %s]',fi,obj.reffeatures(j).features(fi).tostring('mztarget',obj.mztarget(ind,k),'timetarget',obj.time(ind,k),'intensitytarget',obj.fsens(j)*obj.tsens(ind,k),'details',false,'fixedwidth',true),label);
    end
    
    function res=isocheck(obj,mzdata,varargin)
    % Check isotope patterns using original MS files
    % rp is resolving power -- used to flag missing peaks where no other peaks are within mz/rp
      defaults=struct('minpeak',1000,'minisoic',600,'trace',[],'mztol',obj.MZFUZZ,'rp',20000,'noiserel',0.2,'noiseabs',300,'outlierfrac',0.9);
      args=processargs(defaults,varargin);

      res=[];
      for i=1:length(obj.names)
        if isnan(obj.meantime(i))
          continue;
        end
        adduct=obj.astats(i).adduct;
        formula=obj.formula{i};
        if strcmp(obj.ADDUCTS(adduct).name,'M+K')
          % Potassium has a significant ion
          formula(end+1)='K';
        end
        isotopes=Chem.getisotopes(formula,'minabundance',args.minisoic/1e7);
        EMASS=0.00054858;
        for m=1:length(isotopes)
          if strcmp(obj.ADDUCTS(adduct).name,'M+K')
            % Already have the potassium in the formula, just remove an electron mass
            isotopes(m).mass=isotopes(m).mass-EMASS;
          else
            isotopes(m).mass=isotopes(m).mass+obj.ADDUCTS(adduct).mass;
          end
        end
        ic=squeeze(obj.ic(i,adduct,:));
        fsel=find(ic>=args.minpeak & obj.contains(i,:)');
        if isempty(fsel)
          continue;
        end
        if ismember(i,args.trace)
          ti=sprintf('isocheck-%s[%s]',obj.names{i},obj.ADDUCTS(adduct).name);
          setfig(ti);clf;
          tl=tiledlayout('flow');
          title(tl,ti);
          linked=[];
        end
        fres=[];
        for p=1:length(fsel)
          j=fsel(p);
          reftime=obj.time(i,adduct,j);
          rt=interp1(obj.maps{j}.time(:,1),obj.maps{j}.time(:,2),obj.time(i,adduct,j),'linear','extrap');
          [~,scan]=min(abs(mzdata{j}.time-rt));
          peaks=mzdata{j}.peaks{scan};
          maxabund=max([isotopes.abundance]);
          assert(abs(isotopes(1).mass == obj.mztarget(i,adduct))<1e-4);  % Verify that we used the max abund monoisotopic mass
          ictotal=ic(j)/maxabund;
          % Only keep isotopes that are expected to have ic >= minisoic
          keep=find([isotopes.abundance]*ictotal >= args.minisoic);
          if length(keep)<2
            % Nothing other than primary expected
            continue;
          end
          isokeep=isotopes(keep);
          offset=0;   % Correct for shifts in monoisotopic mass
          for k=1:length(isokeep)
            [pkdist,ind]=min(abs(peaks(:,1)-isokeep(k).mass-offset));
            if pkdist<args.mztol
              pval=peaks(ind,2);
            else
              expect=ictotal*isokeep(k).abundance;
              if pkdist>isokeep(k).mass/args.rp && expect>1e4
                fprintf('Missing peak at compound %d, samp %d, iso %d: expected %.0f, closest peak (%.0f=%.0fx) is %.3f away\n', i, j, k, expect,peaks(ind,2),peaks(ind,2)/expect, peaks(ind,1)-isokeep(k).mass-offset);
              end
              pval=0;
            end
            isokeep(k).obsmass=peaks(ind,1);
            isokeep(k).obsic=peaks(ind,2);
            isokeep(k).obs=pval/ictotal;
            isokeep(k).outlier=pval<(isokeep(k).abundance*ictotal*(1-args.noiserel)-args.noiseabs);
            if k==1
              offset=isokeep(k).obsmass-isokeep(k).mass;
            end
            if k==3 && isokeep(k).obs/isokeep(k).abundance<0.8 && ictotal*isokeep(k).abundance>1e5 && pval>0
              fprintf('Low abundance isotope 3 at compound %d, samp %d\n', i,j);
            end
          end
          fres=[fres,struct('sample',j,'ictotal',ictotal,'isotopes',isokeep,'offset',offset,'noutlier',sum([isokeep.outlier]))];
          
          if ismember(i,args.trace)
            fprintf('%s[%s] in %s, rt=%.2f, ic=%.0f\n', obj.names{i}, obj.ADDUCTS(adduct).name,obj.samples{j},rt,ic(j));
            nexttile;
            semilogy(peaks(:,1),peaks(:,2),'o-');
            hold on;
            plot([isokeep.mass],[isokeep.abundance]*ic(j)/isokeep(1).abundance,'x');
            ax=axis; ax(1:2)=[min([isokeep.mass])-1,max([isokeep.mass])+1]; ax(4)=ic(j)*1.2; axis(ax);
            xlabel('m/z');
            ylabel('ion count');
            title(sprintf('%s - %s',obj.samples{j},isokeep(1).name));
            linked(end+1)=gca;
          end
        end
        if ~isempty(fres)
          res=[res,struct('compound',i,'adduct',adduct,'samples',fres,'noutlier',sum([fres.noutlier]>0))];
        end
        if ismember(i,args.trace)
          linkaxes(linked,'x');
        end
      end

      % List outliers
     outlierfrac=arrayfun(@(z) z.noutlier/length(z.samples), res);
      if any(outlierfrac>=args.outlierfrac)
        fprintf('Outliers in at least %.0f%% of files:\n',args.outlierfrac*100);
        for i=1:length(res)
          if outlierfrac(i)>=args.outlierfrac
            r=res(i);
            fprintf('%4d %s[%s] %s %.5f noutliers=%d\n', r.compound, obj.names{r.compound}, obj.ADDUCTS(r.adduct).name,r.samples(1).isotopes(1).name, r.samples(1).isotopes(1).mass, r.noutlier);
            for j=1:length(r.samples)
              s=r.samples(j);
              fprintf(' %2d %-10s %6.0f', s.sample, obj.samples{s.sample},s.ictotal)
              for k=2:length(s.isotopes)
                iso=s.isotopes(k);
                fprintf(' %5.0f/%5.0f',[iso.obs,iso.abundance]*s.ictotal);
                if iso.outlier
                  fprintf('*');
                else
                  fprintf(' ');
                end
              end
              fprintf('\n');
            end
          end
        end
      else
        fprintf('No compounds with outliers in at least %.0f%% of the files\n',args.outlierfrac*100);
      end
      
      % Plot observed vs. expected ion counts and outlier dependence
      setfig('isocheck');clf;
      data=[];
      for i=1:length(res)
        for k=1:length(res(i).samples)
          r=res(i).samples(k);
          for j=2:length(r.isotopes)
            iso=r.isotopes(j);
            data(end+1,:)=r.ictotal*[iso.abundance,iso.obs];
          end
        end
      end
      loglog(data(:,1),data(:,2),'.','MarkerSize',1);
      thresh=data(:,1)*(1-args.noiserel)-args.noiseabs;   % Points below this are outliers
      thresh(thresh<=0)=nan;
      hold on;
      outlier=data(:,2)<thresh;
      loglog(data(outlier,1),data(outlier,2),'r.','MarkerSize',1);
      
      % Plot threshold
      hold on;
      ax=axis;
      [~,ord]=sort(data(:,1));
      plot(data(ord,1),thresh(ord),'r');
      xlabel('Expected ion count');
      ylabel('Observed ion count');
      title(sprintf('isocheck noise=%f*ic+%f',args.noiserel,args.noiseabs));
      
      ax=axis; ax(3:4)=[100,1e6];axis(ax);

      % Plot observed vs. expected ion counts for 3rd isotope of compounds with Sulfur
      setfig('isocheck-Sulfur');clf;
      tiledlayout('flow');
      data=[];
      for i=1:length(res)
        hassulfur=any(res(i).samples(1).isotopes(1).name(end)=='S');
        for k=1:length(res(i).samples)
          r=res(i).samples(k);
          if length(r.isotopes)>=3
            iso=r.isotopes(3);
            data(end+1,:)=[r.ictotal*[iso.abundance,iso.obs],hassulfur];
          end
        end
      end
      has=data(:,3)~=0;
      nexttile;
      loglog(data(has,1),data(has,2),'.','MarkerSize',1);
      hold on;
      ax=axis; 
      plot(ax(1:2),ax(1:2),'r:');
      xlabel('Expected ion count');
      ylabel('Observed ion count');
      title(sprintf('isocheck Sulfur#3 only'));
      nexttile;
      sel=data(:,1)>2e4 & has;
      pdfplot(data(sel,2)./data(sel,1));
      hold on;
      sel=data(:,1)>2e4 & ~has;
      pdfplot(data(sel,2)./data(sel,1));
      legend('Has Sulfur','No Sulfur');
      xlabel('Obs/Expected');
      set(gca,'XScale','log');
      ax=axis;ax(2)=2;axis(ax);
      
      % Plot mass differential
      setfig('isocheck-mass');clf;
      tiledlayout('flow');
      data=[];
      for i=1:length(res)
        for j=1:length(res(i).samples)
          r=res(i).samples(j);
          iso1=r.isotopes(1);
          for k=2:length(r.isotopes)
            iso=r.isotopes(k);
            data(end+1,:)=[k,iso.obsic,iso.mass,iso.obsmass,iso1.obsmass-iso1.mass];
          end
        end
      end
      leg={};
      nexttile;
      for k=2:max(data(:,1))
        sel=data(:,1)==k;
        semilogx(data(sel,2),data(sel,4)-data(sel,3)-data(sel,5),'.','MarkerSize',1);
        leg{end+1}=sprintf('Isotope %d',k);
        hold on;
      end
      set(gca,'XScale','log');
      ax=axis;
      semilogx(ax(1:2),args.mztol*[1,1],'r:');
      semilogx(ax(1:2),-args.mztol*[1,1],'r:');
      ax(3:4)=args.mztol*[-1,1]*1.5;
      axis(ax);
      xlabel('Ion Count');
      ylabel('Mass Error');
      legend(leg);
      title('Mass Error');

      nexttile;
      leg={};
      dmz=data(:,4)-data(:,3)-data(:,5);
      for k=2:max(data(:,1))
        sel=data(:,1)==k & data(:,2)>10000 & abs(dmz)<=args.mztol;
        if sum(sel)>100
          pdfplot(dmz(sel));
          hold on;
          leg{end+1}=sprintf('Isotope %d',k);
        end
      end
      xlabel('Mass Error');
      legend(leg);
      
    end
    
    function getinfo(obj,name,varargin)
      defaults=struct('mzdata',[],'adduct',[],'falsethresh',[]);
      args=processargs(defaults,varargin);

      if ischar(name)
        ind = obj.find(name);
      else
        ind=name;
      end
      if isempty(args.adduct)
        % Check which adduct was used to label, if any
        if ~isempty(obj.astats) && ind<=length(obj.astats) && ~isempty(obj.astats(ind).adduct)
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
      fprintf('%s[%s] (%d): m/z=%8.4f t=%.2f [%.2f-%.2f] sens=%.0f',obj.names{ind},obj.ADDUCTS(k).name, ind, obj.mztarget(ind,k),obj.meantime(ind),obj.timewindow(ind,:),obj.tsens(ind,k));
      if ~isempty(obj.astats) && ind<=length(obj.astats) && ~isempty(obj.astats(ind).run) && obj.astats(ind).adduct==args.adduct
        s=obj.astats(ind);
        fprintf(', FP=%d, FN=%d, hits=%dg,%dL,%dH, miss=%ds,%dw',s.FP,s.FN,s.hitgood,s.hitlow,s.hithigh,s.missstrong,s.missweak);
      end
      fprintf('\n');
      if ~isempty(args.mzdata)
        setfig([obj.names{ind},'-',obj.ADDUCTS(args.adduct).name]);
        t=tiledlayout('flow');
        title(t,sprintf('%s m/z=%.4f t=%.2f',obj.names{ind},obj.mztarget(ind,k),obj.meantime(ind)));
      end
      
      for j=1:length(obj.files)
        if ~obj.contains(ind,j)
          continue;
        end
        fprintf('%2d %-15.15s %5.2f',j,obj.samples{j},obj.fsens(j));

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
        if isempty(args.falsethresh)
          % Use the same value that was used during assignTimes
          args.falsethresh=obj.astats(ind).args.falsethresh;
        end
        firstfalse=true;
        for j=1:length(obj.files)
          fi=obj.featureindex(ind,k,j);
          if isfinite(fi) && ~obj.contains(ind,j)
            if isnan(obj.normic(ind,k,j))
              fprintf('Normalized IC not set, run checksensitivity to get false positives\n');
              firstfalse=false;
              break;
            end
            if obj.normic(ind,k,j)>=args.falsethresh
              if firstfalse
                fprintf('False positives with normIC >= %.3f:\n',args.falsethresh);
                firstfalse=false;
              end
              fprintf(' %2d %-14.14s %5.2f%s\n',j,obj.samples{j},obj.fsens(j),obj.featstring(ind,j,k,fi));
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
        if firstfalse
          fprintf('No false positives with normIC >= %.3f\n', args.falsethresh);
        end
      end
      % Aliases
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
            fprintf('\t%-6.6s[%-4.4s]: m/z=%8.4f (d=%4.0f) t=%.2f (d=%.2f) meanic=%.0f, in files [%s]\n',...
                    obj.names{i},obj.ADDUCTS(k1).name, obj.mass(i)+obj.ADDUCTS(k1).mass,...
                    (obj.mass(i)+obj.ADDUCTS(k1).mass-(obj.mass(ind)+obj.ADDUCTS(k).mass))*1e4,...
                    obj.meantime(i),obj.meantime(i)-obj.meantime(ind),imeanic,...
                    strjoin(arrayfun(@(z) sprintf('%d',z), find(obj.contains(i,:)),'Unif',false),','));
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
        set(gca,'XTick',round(obj.timewindow(ind,:),1));
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
          k1=obj.astats(ii).adduct;
          m1=(obj.mass(ii)+obj.ADDUCTS(k1).mass);
          k2=obj.astats(jj).adduct;
          m2=(obj.mass(jj)+obj.ADDUCTS(k2).mass);
          mdiff=m1-m2;
          if abs(mdiff) < 2*obj.MZFUZZ
            fprintf('Overlap: %6.6s+%-4.4s (%.4f,[%.2f,%.2f]) and %6.6s+%-4.4s (%.4f,[%.2f,%.2f]) dmz=%2.0f dt=%.2f\n', obj.names{ii},obj.ADDUCTS(k1).name,m1,obj.timewindow(ii,:),obj.names{jj},obj.ADDUCTS(k2).name,m2,obj.timewindow(jj,:),abs(mdiff)*1e5,abs(obj.meantime(ii)-obj.meantime(jj)));
            nover=nover+1;
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
    
    function export(obj,filename,varargin)
    % Export via CSV (can be used as import to mzmine)
    % Each row has fields:  ID, m/z, retention time, identity, formula, tsens
      defaults=struct('adducts',[],'notimes',false,'ind',[]);
      args=processargs(defaults,varargin);
    
      if isempty(args.ind)
        args.ind=true(size(obj.mass));
      end
      if islogical(args.ind)
        args.ind=find(args.ind);
      end
      if isempty(args.adducts)
        args.adducts=1:length(obj.ADDUCTS);
      end
      fd=fopen(filename,'w');
      fprintf(fd,'ID, m/z, retention time, identity, formula, tsens\n');
      for ii=1:length(args.ind)
        i=args.ind(ii);
        if isempty(args.adducts)
          % Export only best one (highest IC)
          bestadduct=1;
          for j=1:length(args.adducts)
            if obj.tsens(i,j)>obj.tsens(i,bestadduct)
              bestadduct=j;
            end
          end
          fprintf(fd,'%d,%.4f,%.2f,%s[%s],%s,%f\n', i, obj.mass(i)+obj.ADDUCTS(bestadduct).mass,t/60.0,obj.names{i},obj.ADDUCTS(bestadduct).name,obj.formula{i},obj.tsens(i,bestadduct));
        else
          for j=1:length(args.adducts)
            t=obj.meantime(i);
            if isnan(t) || args.notimes
              t=0;
            end
            fprintf(fd,'%d,%.4f,%.2f,%s[%s],%s\n', i, obj.mass(i)+obj.ADDUCTS(j).mass,t/60.0,obj.names{i},obj.ADDUCTS(j).name,obj.formula{i});
          end
        end
      end
      fclose(fd);
    end
    
  end
end
