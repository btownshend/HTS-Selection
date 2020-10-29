classdef MassSpec < handle
  properties
    path;    % Full path
    name;    % Short name for this file
    moles;   % Nominal number of moles (per compound) loaded
    peaks;  % Peaks from mzxml2peaks;  peaks{i}(k,1:2) is [mz,ic] for elution i, at peak k 
    time;   % Time of elutions (in minutes)
    mzrange;   % Range [low,high] of m/z to plot/consider
    featurelists;   % Array of feature lists
    mzxml;
  end
  
  properties(Transient)
    resamp;   % struct (mz,y,n) of uniformly sampled data
  end
  
  properties(Constant)
    mzres=0.01;
  end
  
  methods(Static)
    function c=comparepeaks(objs,varargin)
    % Compare peak lists
      defaults=struct('mztol',0.01,'timerange',[],'minic',5000,'mz',[],'nlist',20,'ref',1,'maxtiles',12);
      args=processargs(defaults,varargin);
      
      c=[];
      for i=1:length(objs)
        fprintf('Locating peaks in %s\n', objs{i}.name);
        pks{i}=objs{i}.allpeaks('mztol',args.mztol,'timerange',args.timerange,'minic',args.minic,'mz',args.mz);
        c(end+1:end+size(pks{i},1),[1,i+1])=pks{i};
      end

      % Sort peak list by descending IC
      [~,ord]=sort(sum(c(:,2:end),2),'desc');
      c=c(ord,:);
      
      % Merge peaks which differ by <mztol
      fprintf('Merging all peaks from %d objs...',length(objs));
      for i=1:size(c,1)
        if isfinite(c(i,1))
          close=find(abs(c(i,1)-c(:,1))<args.mztol);
          c(i,2:end)=sum(c(close,2:end),1);
          c(close(close~=i),1)=nan;
        end
      end
      fprintf('done\n');
      % Remove merged peaks
      c=c(isfinite(c(:,1)),:);
      
      ntiles=0;
      for i=1:length(objs)
        if i==args.ref
          continue;
        end
        if mod(ntiles,args.maxtiles)==0
          setfig(sprintf('comparepeaks:%d',i));clf;
          tiledlayout('flow');
        end
        ntiles=ntiles+1;
        nexttile;
        ratio=c(:,args.ref+1)./c(:,i+1);
        ratio(isinf(ratio))=nan;
        ti=sprintf('Ratio Med:%.2f Mean:%.2f CV=%.2f\n', nanmedian(ratio),nanmean(ratio),nanstd(ratio)./nanmean(ratio));
        loglog(c(:,args.ref+1),c(:,i+1),'.');
        hold on;
        rng=[args.minic,max(max(c(:,2:end)))];
        plot(rng,rng,':r');
        xlabel(objs{args.ref}.name,'Interpreter','none');
        ylabel(objs{i}.name,'Interpreter','none');
        title(ti,'interpreter','none');
        fprintf('%s: %s\n', objs{i}.name, ti);
      end
      for i=1:min(args.nlist,size(c,1))
        fprintf('%8.4f %s\n',c(i,1),sprintf('%4.2f ',c(i,2:end)/c(i,args.ref+1)));
      end
    end
  end
  
  methods
    function obj=MassSpec(path,varargin)
      defaults=struct('debug',false,'savexml',false);
      args=processargs(defaults,varargin);
      fprintf('Loading %s...\n', path);
      mzxml=mzxmlread(path);
      [obj.peaks,obj.time]=mzxml2peaks(mzxml);
      good=true(size(obj.peaks));
      for i=1:length(obj.peaks)
        maxp=max(obj.peaks{i}(:,2));
        if maxp>1e10
          fprintf('Bad scan %d:  max peak is %g - deleting\n', i, maxp);
          good(i)=false;
        end
      end
      if any(~good)
        obj.peaks=obj.peaks(good);
        obj.time=obj.time(good);
      end
      obj.time=obj.time/60;   % Convert to minutes
      obj.mzrange=[min(cellfun(@(z) z(1,1), obj.peaks)),max(cellfun(@(z) z(end,1), obj.peaks))];
      obj.path=path;
      z=strsplit(obj.path,'/');
      obj.name=z{end};
      if args.savexml
        obj.mzxml=mzxml;
      end
    end

    function setLoad(obj,moles)
      obj.moles=moles;
    end
    
    function filter(obj,trange,mzrange)
      if nargin>=2 && ~isempty(trange)
        sel=find(obj.time>=trange(1) & obj.time<=trange(2));
        obj.time=obj.time(sel);
        obj.peaks=obj.peaks(sel);
      end
      if nargin>=3 && ~isempty(mzrange)
        obj.mzrange=mzrange;
      end
    end

    function keepmz(obj,mz,varargin)
    % Keep only the peaks within mztol of an entry in mz
      defaults=struct('mztol',0.01);
      args=processargs(defaults,varargin);
      
      fprintf('Keeping only peaks within %f of one of %d m/z values...',args.mztol, length(mz));
      nfinal=0; norig=0;
      mz=reshape(mz,1,[]);
      for i=1:length(obj.peaks)
        if i==1
          fprintf('%d/%d...',i,length(obj.peaks));
        elseif mod(i,100)==0
          fprintf('%d...',i);
        end
          
        
        p=obj.peaks{i};
        norig=norig+size(p,1);
        tic;sel=ismembertol(p(:,1),mz,args.mztol,'DataScale',1);
        if toc>1
          keyboard;
        end
        nfinal=nfinal+sum(sel);
        obj.peaks{i}=p(sel,:);
      end
      fprintf('done, retained %d/%d peaks (%.1f%%)\n', nfinal, norig,nfinal/norig*100);
    end
    
    
    function resample(obj)
      n=round(diff(obj.mzrange)/obj.mzres+1);
      fprintf('Resampling peaks...');
      [mz,y]=msppresample(obj.peaks,n,'Range',obj.mzrange,'FWHH',obj.mzres);
      fprintf('done\n');
      obj.resamp=struct('mz',mz,'y',y,'n',n);
    end

    function heatmap(obj,markers)
      if isempty(obj.resamp)
        obj.resample();
      end
      if nargin<2
        markers=[];
      end
      if isempty(markers)
        msheatmap(obj.resamp.mz,obj.time,log(obj.resamp.y),'resolution',obj.mzres);
      else
        msheatmap(obj.resamp.mz,obj.time,log(obj.resamp.y),'resolution',obj.mzres,'markers',markers);
      end
    end

    function [ic,mz,time]=mzscan(obj, mztarget, varargin)
    % Get ion counts for given mztarget across elution times
    % Also return weighted mean m/z for each elution time
      defaults=struct('mztol',0.01,'timerange',[-inf,inf]);
      args=processargs(defaults,varargin);
      % Build a table of elution vs IC for this M/Z
      first=find(obj.time>=args.timerange(1),1);
      last=find(obj.time<=args.timerange(2),1,'last');

      ic=zeros(last-first+1,1);
      mz=nan(last-first+1,1);
      time=obj.time;time=time(first:last);
      for i=first:last
        p=obj.peaks{i};
        pks=p(:,1);
        sel=abs(pks-mztarget)<=args.mztol;
        if any(sel)
          ic(i-first+1)=sum(p(sel,2));
          mz(i-first+1)=sum((p(sel,1).*p(sel,2)))/ic(i-first+1);
        end
      end
    end
    
    function pks=allpeaks(obj,varargin)
    % Get peaks list summed across all elution times and merged with mztol
      defaults=struct('mztol',0.01,'timerange',[],'minic',1000,'mz',[]);
      args=processargs(defaults,varargin);
      % Build a table of elution vs IC for this M/Z
      if isempty(args.timerange)
        first=1;
        last=length(obj.peaks);
      else
        first=find(obj.time>=args.timerange(1),1);
        last=find(obj.time<=args.timerange(2),1,'last');
      end
      % Merge all peaks
      fprintf('Concatenating peaks across %d elutions\n', last-first+1);
      pks=[];
      if isempty(args.mz)
        for i=first:last
          p=obj.peaks{i};
          pks=[pks;p(p(:,2)>args.minic/10,:)];
        end
      else
        for i=first:last
          p=obj.peaks{i};
          sel=min(abs(p(:,1)-args.mz(:)'),[],2)<=args.mztol;
          pks=[pks;p(sel & p(:,2)>args.minic/10,:)];
        end
      end
      fprintf('Have a total of %d peaks\n', size(pks,1));

      if ~isempty(args.mz)
        % Set them to specified mz
        for i=1:size(pks,1)
          [~,ind]=min(abs(args.mz-pks(i,1)));
          pks(i,1)=args.mz(ind);
        end
      end

      % Sort peak list by descending IC
      [~,ord]=sort(pks(:,2),'desc');
      pks=pks(ord,:);
      % Merge peaks which differ by <mztol
      fprintf('Merging...');
      todo=1:length(pks);
      mpks=[];
      i=1;
      while ~isempty(todo)
        if mod(i,1000)==0
          fprintf('%.0f%%...',(1-length(todo)/size(pks,1))*100);
        end
        close=abs(pks(todo(1),1)-pks(todo,1))<args.mztol;
        mpks(end+1,:)=[pks(todo(1),1),sum(pks(todo(close),2))];
        todo=todo(~close);
        i=i+1;
      end
      fprintf('done\n');
      % Remove merged peaks
      pks=mpks;
      fprintf('Have %d merged peaks\n', size(pks,1));

      % Keep only ones with >= minic
      pks=pks(pks(:,2)>=args.minic,:);
      fprintf('Have %d merged peaks with IC >= %.0f\n', size(pks,1),args.minic);
    end
      
    function res=findcompound(obj, mztarget, varargin)
    % Find compound directly (not using chromatograms)
    % Returns struct containing possible peaks (integrated over peak in both time and m/z) at different elutions times
    % Only peaks with ion count >= peakratio* maximum peak ion count are returned
    % If sdf provided, only used to construct a name for the compound
    % overseg is minimum separation between peaks; could be less than TIMEFUZZ, since the fuzz also includes misalignment of times
      defaults=struct('elutetime',[],'mztol',0.01,'timetol',30/60,'debug',false,'ignoreelutetimes',[],'peakratio',0.05,'overseg',nan,'showplot',false,'heightfilter',100);
      args=processargs(defaults,varargin);
      
      if isnan(args.overseg)
        args.overseg=args.timetol;
      end
      
      % Build a table of elution vs IC for this M/Z
      if isempty(args.elutetime)
        [ic,mz,time]=obj.mzscan(mztarget,'mztol',args.mztol);
      else
        [ic,mz,time]=obj.mzscan(mztarget,'mztol',args.mztol,'timerange',[args.elutetime-args.timetol,args.elutetime+args.timetol]);
      end

      % Handle any other time limitations given on command line
      if ~isempty(args.ignoreelutetimes)
        sel=~(min(abs(time-args.ignoreelutetimes(:)'),[],2)<args.timetol);
        time=time(sel);
        ic=ic(sel);
        mz=mz(sel);
      end

      
      % Find peaks
      maxtime=[]; maxic=[]; maxmz=[]; pfwhh=[];
      if ~isempty(ic) && length(ic)>1 && any(ic>0)
        % mspeaks won't find peaks that are only in the last sample (or probably first)
        [p,pfwhh]=mspeaks([time(1)-1;time;time(end)+1],[0;ic;0],'denoising',false,'OVERSEGMENTATIONFILTER',args.overseg,'heightfilter',args.heightfilter,'showplot',args.showplot,'style','fwhhtriangle');
        if ~isempty(p)
          for i=1:size(p,1)
            sel=time>=pfwhh(i,1) & time<=pfwhh(i,2);
            maxic(i)=sum(ic(sel));   % Integrate instead of just using the max height
            maxmz(i)=nansum(mz(sel).*ic(sel))/maxic(i);
            maxtime(i)=p(i,1);
          end
        else
          fprintf('m/z scan for %f gave %d peaks, but none from mspeaks\n', mztarget, sum(ic>0));
          sel=ic>=args.heightfilter;
          maxmz=mz(sel)';
          maxtime=time(sel)';
          maxic=ic(sel)';
          for i=1:length(time)
            pfwhh(i,1:2)=[1,1]*time(i);
          end
        end
        [maxic,ord]=sort(maxic,'desc');
        maxmz=maxmz(ord);
        maxtime=maxtime(ord);
        pfwhh=pfwhh(ord,:);
        keep=maxic>=max(maxic)*args.peakratio;
        maxic=maxic(keep);
        maxmz=maxmz(keep);
        maxtime=maxtime(keep);
        pfwhh=pfwhh(keep,:);
      end
      assert(all(isfinite(maxmz)));
      desc=sprintf('m/z=%.4f+=%.3f',mztarget,args.mztol);
      if ~isempty(args.elutetime)
        desc=[desc,sprintf(' T=%.1f+=%.1f',args.elutetime,args.timetol)];
      elseif ~isempty(args.ignoreelutetimes)
        desc=[desc,sprintf(' T!=[%s]+=%.1f',sprintf('%.1f,',args.ignoreelutetimes),args.timetol)];
      end

      res=struct('mztarget',mztarget,'desc',desc,'findargs',args,'mz',maxmz,'time',maxtime,'ic',maxic,'pfwhh',pfwhh);
      if isempty(maxic)
        if args.debug
          fprintf('No match to %s\n',desc);
        end
      elseif length(maxmz)==1
        if args.debug
          fprintf('Unique peak for %s:  T=%.2f [%.2f-%.2f] IC=%.0f\n',desc, maxtime(1),pfwhh(1,:),maxic(1));
        end
      else
        if args.debug
          fprintf('Have %d peaks for %s with peak ratio=%f\n',length(maxmz), desc, maxic(2)/maxic(1));
        end
      end
    end

    function buildchromatograms(obj,varargin)
    % Build chromatograms using ADAP algorithm (as documented in mzmine2 help)
    % We define ε to be the mass tolerance parameter m/z.
    % 1. Take all the data points in a data file, sort them by their intensities, and remove those points (mostly noise) below a certain 
    %    intensity threshold.
    % 2. Starting with the most intense data point, the first EIC is created.
    % 3. For this EIC, establish an immutable m/z range that is the data point's m/z plus or minus ε where ε is specified by the user.
    % 4. The next data point, which will be the next most intense, is added to an existing EIC if its m/z value falls within its m/z range.
    % 5. If the next data point does not fall within an EICs m/z range, a new EIC is created. New EICs are only created if the point meets 
    %    the minimum start intensity requirement set by the user.
    % 6. An m/z range for a new EIC is created the same way as in step (3) except the boundaries will be adjusted to avoid overlapping with 
    %    pre-existing EICs. As an example consider an existing EIC with m/z range (100.000,100.020) for ε=0.01. If the new EIC is initialized 
    %    with a data point having an m/z value of 100.025, then this new EIC will have a m/z range set to be (100.020,100.035) rather than
    %    (100.015, 100.035).
    % 7. Repeat steps (4)-(6) until all the data has been processed.
    % 8. Finally, a post processing step is implemented. Only EICs with a user defined number of continuous points (mingroupsize) above 
    % a user defined intensity threshold (groupthresh) are kept.
      defaults=struct('mztol',0.01,'debug',false,'mingroupsize',2,'groupthresh',500,'minintensity',1000,'noise',500);
      args=processargs(defaults,varargin);
      % Build list of all peaks with time index as 3rd column
      allpks=[];
      for i=1:length(obj.peaks)
        p=obj.peaks{i};
        p=p(p(:,2)>=args.noise & p(:,1)>=obj.mzrange(1) & p(:,1)<=obj.mzrange(2),:);
        p(:,3)=i;
        allpks=vertcat(allpks,p);
      end
      if args.debug
        fprintf('Have a total of %d peaks with IC >= %.0f\n', size(allpks,1), args.noise);
      end
      
      % Sort by level
      [~,ord]=sort(allpks(:,2),'descend');
      allpks=allpks(ord,:);
      % Buld EIC
      eicranges=zeros(0,2);   % [mzlow, mzhigh]
      fl=FeatureList([obj.name,' chromatograms'],'buildchromatogram',args);
      for i=1:size(allpks,1)
        if allpks(i,2)>=args.minintensity
          mzrange=[allpks(i,1)-args.mztol,allpks(i,1)+args.mztol];
          % Reduce mass range if it overlaps
          mzrange(1)=max([eicranges(eicranges(:,2)<allpks(i,1),2);mzrange(1)]);
          mzrange(2)=min([eicranges(eicranges(:,2)>allpks(i,1),2);mzrange(2)]);
          % Pull out all other peaks that will match
          sel=find(allpks(i:end,1)>=mzrange(1) & allpks(i:end,1)<=mzrange(2) & allpks(i:end,2)>0)+i-1;
          % Check if meets the the group thresholds
          gscans=allpks(sel(allpks(sel,2)>=args.groupthresh),3);
          gscans=sort(unique(gscans));
          ok=any(gscans(args.mingroupsize:end)-gscans(1:end-args.mingroupsize+1)==args.mingroupsize-1);
          if args.debug
            ngroup=0;
            for k=1:length(gscans)
              nn=1;
              for m=k+1:length(gscans)
                if gscans(m)-gscans(k) == m-k
                  nn=nn+1;
                else
                  break;
                end
              end
              if nn>ngroup
                ngroup=nn;
              end
            end
            fprintf('m/z=%.4f, IC=%.0f new EIC [%.4f,%.4f] with %d scans (%d>%.0f, group %d) pass=%d\n', ...
                    allpks(i,1:2), mzrange, length(sel),length(gscans),args.groupthresh,ngroup, ok);
          end
          if (ok)
            % New EIC feature
            eicranges(end+1,:)=mzrange;
            p=nan(length(obj.time),3);
            ap=allpks(sel,:);
            for j=1:length(obj.time)
              ind=ap(:,3)==j;
              if any(ind)
                p(j,:)=[mean(ap(ind,1)),sum(ap(ind,2)),j];
              else
                p(j,:)=[nan,0,j];
              end
            end
            p(:,3)=obj.time(p(:,3)); % Convert to time (but not earlier since we use for gscans)
            area=sum(p(:,2));
            feature=Feature(p,sprintf('%.4f-%.4f',mzrange));
            fl.append(feature);
          end
          % Blank out the consumed ones (whether above group threshold or not)
          allpks(sel,2)=0;
        end
      end
      if args.debug
        fprintf('Found %d distinct EICS\n', size(eicranges,1));
      end
      % Sort by mz
      fl.sortbymz();
      % Append
      obj.featurelists=[obj.featurelists,fl];
    end

    function fl=targetedFeatureDetect(obj,mz,varargin)
    % Build chromatograms using targets
    % Retention time arg is optional
      defaults=struct('mztol',0.01,'timetol',[],'debug',false,'rt',[],'names',{{}},'noise',500);
      args=processargs(defaults,varargin);
      
      [mz,ord]=sort(mz);

      if ~isempty(args.rt)
        assert(length(args.rt)==length(mz));
        args.rt=args.rt(ord);
      end
      if ~isempty(args.names)
        assert(length(args.names)==length(mz));
        args.names=args.names(ord);
      end
      mzrange(:,1)=mz-args.mztol;
      mzrange(:,2)=mz+args.mztol;
      % Adjust overlaps
      if isempty(args.rt)
        % Merge duplicate or overlapping m/z
        ndup=0;
        while true
          dup=find(mzrange(1:end-1,2)>mzrange(2:end,1),1);
          if isempty(dup)
            break;
          end
          mzrange(dup,2)=mzrange(dup+1,2);
          mzrange=mzrange([1:dup,dup+2:end],:);
          if ~isempty(args.names)
            args.names{dup}=strjoin(args.names([dup,dup+1]),',');
            args.names=args.names([1:dup,dup+2:end]);
          end
          ndup=ndup+1;
        end
        if ndup>0
          fprintf('Merged %d duplicate m/z\n', ndup);
        end
      else
        % With RT given, need to handle overlaps -- TODO
        for i=1:length(mz)
          for j=i+1:length(mz)
            if mz(j)-mz(i) < 2*args.mztol
              break;   % Distinct m/z
            end
            if abs(args.rt(j)-args.rt(i)) < 2*args.timetol
              fprintf('Overlap: %.4f@%.2f and %.4f@%.2f\n', mz(i), args.rt(i), mz(j), args.rt(j));
            end
          end
        end
        
      end
      
      % Buld EIC for each mass
      fl=FeatureList([obj.name,' chromatograms'],'targetedFeatureDetect',args);
      if args.debug
        fprintf('Extracting features for %d masses...',size(mzrange,1));
      end
      
      for i=1:size(mzrange,1)
        if args.debug && mod(i,500)==1
          fprintf('%d...',i);
        end
        
        if ~isempty(args.rt)
          timerange=args.rt(i)+[-1,1]*args.timetol;
        end
        
        if ~isempty(args.names)
          fname=args.names{i};
        elseif isempty(args.rt)
          fname=sprintf('%.4f-%.4f',mzrange(i,:));
        else
          fname=sprintf('%.4f-%.4f@%.1f-%.1f',mzrange(i,:),timerange);
        end
          
        % New EIC feature
        p=[];
        for j=1:length(obj.time)
          if ~isempty(args.rt) && (obj.time(j)<timerange(1) || obj.time(j)>timerange(2))
            continue;
          end
          peaks=obj.peaks{j};
          ind=peaks(:,1)>=mzrange(i,1) & peaks(:,1)<mzrange(i,2);
          if any(ind)
            psel=peaks(ind,:);
            [~,ia]=max(psel(:,2));
            p(end+1,:)=[psel(ia),sum(psel(:,2)),j];
          else
            p(end+1,:)=[mean(mzrange(i,:)),0,j];
          end
        end
        if isempty(p)
          continue;   % Not found
        end
        if ~any(p(:,2)>args.noise)
          continue;   % Too small
        end
        %p=p(1:find(isfinite(p(:,1)),1,'last'),:);   % Removing trailing nans
        first=max(1,find(p(:,2)>0,1)-1);
        last=min(size(p,1),find(p(:,2)>0,1,'last')+1);
        p=p(first:last,:);
        p(:,3)=obj.time(p(:,3)); % Convert from scan to time
        feature=Feature(p,fname);
        fl.append(feature);
      end
      if args.debug
        fprintf('done\n');
      end
      % Sort by mz
      fl.sortbymz();
      if nargout==0
        % Append
        obj.featurelists=[obj.featurelists,fl];
      end
    end

    function deconvolve(obj,varargin)
      defaults=struct('src',[]);
      args=processargs(defaults,varargin);

      if isempty(args.src)
        args.src=find(~strcmp({obj.featurelists.src},'deconvolve'),1,'last');
      end
      fprintf('Deconvolve of %s\n', obj.featurelists(args.src).name);
      assert(isempty(strfind(obj.featurelists(args.src).name,'deconvoluted')));  % No double deconvolutes
      fl=obj.featurelists(args.src).deconvolve(varargin{:});
      obj.featurelists=[obj.featurelists,fl];
    end      

    function listfeatures(obj)
      for i=1:length(obj.featurelists)
        f=obj.featurelists(i);
        fprintf('%d name:%-40s  src:%-17s features:%d\n', i, f.name, f.src, length(f.features));
      end
    end
    
    function ploteic(obj,mz,varargin)
      defaults=struct('debug',false,'mztol',0.01,'newfig',true,'timerange',[-inf,inf]);
      args=processargs(defaults,varargin);

      if isempty(obj.featurelists)
        fprintf('Warning - no chromatograms; need to buildchromatograms()\n');
      elseif length(obj.featurelists)<2
        fprintf('Warning - no deconvolved chromatograms\n');
      end
      
      ti=sprintf('%s %.4f - %.4f',obj.name, mz-args.mztol,mz+args.mztol);
      if args.newfig
        setfig(ti);clf;
      end

      % Plot base traces from raw peaks
      [ic,~,time]=obj.mzscan(mz,'mztol',args.mztol,'timerange',args.timerange);
      plot(time,ic);
      legend('Base peaks','Location','best');
      hold on;
      
      %sel=find(arrayfun(@(z) z.mzrange(1)-args.mztol<=mz && z.mzrange(2)+args.mztol>=mz,obj.eic));
      for ifl=1:length(obj.featurelists)
        fl=obj.featurelists(ifl);
        fl.ploteic(mz,'mztol',args.mztol,'prefix',sprintf('%d',ifl),'timerange',args.timerange);
      end
      xlabel('Time');
      ylabel('Intensity');
      title(ti);
    end
    
    function plotmz(obj,time,varargin)
      defaults=struct('maxlabels',10,'newfig',true);
      args=processargs(defaults,varargin);

    % Plot m/z vs ioncount for trace nearest given time
      [~,closest]=min(abs(time-obj.time));
      pks=obj.peaks{closest};
      ti=sprintf('m/z @ T=%.2f (%d)',obj.time(closest),closest);
      if args.newfig
        setfig(ti);clf;
      end
      stem(pks(:,1),pks(:,2),'.');
      title(ti)
      xlabel('m/z')
      ylabel('Relative Intensity')
      hold on;
      ax=axis;
      % List major peaks
      [~,ord]=sort(pks(:,2),'desc');
      for ii=1:args.maxlabels
        i=ord(ii);
        fprintf('M/Z=%.4f, IC=%6.0f\n', pks(i,:));
        h=text(pks(i,1),pks(i,2),sprintf('%.4f',pks(i,1)));
        set(h,'VerticalAlignment','bottom','HorizontalAlignment','left');
      end
    end
    
    function ic=TIC(obj)
      ic=cellfun(@(z) sum(z(:,2)),obj.peaks);
    end
    
    function plotTIC(obj,newfig,time)
    % Plot TIC showing best traces
      ti=[obj.name,'-TIC'];
      if nargin<2 || isempty(newfig) || newfig
        setfig(ti);clf;
      end
      if nargin<3 || isempty(time)
        time=obj.time;
      end
      plot(time,obj.TIC());
      title(ti,'Interpreter','None')
      xlabel('Retention Time (min)')
      ylabel('Relative Intensity')
      hold on;
      ax=axis;
    end

    function map=timealign(obj,o2,step,window)
    % Find mapping between elution times
      t1=obj.TIC();
      t2=o2.TIC();
      if nargin<3
        step=10;
      end
      if nargin<4
        window=500;
      end
      map=[];
      for i=1:step:length(t1)-window
        ind=finddelay(t1(i:i+window-1).*hamming(window),t2(i:i+window-1).*hamming(window),window/10);
        map(end+1,:)=[i+window/2,i+window/2+ind];
      end
      setfig(sprintf('timealign %s-%s',obj.name,o2.name));
      plot(map(:,1),map(:,2)-map(:,1),'.-');
      xlabel(obj.name);
      ylabel(o2.name);
    end
        
    function findtarget(obj,formula,varargin)
    % Scan for a particular m/z include adduct and isotope search
      defaults=struct('dbsave',false,'mztol',0.0005,'timetol',0.6,'noise',500,'debug',false,...
                      'adducts',struct('name',{'M+H','M+Na','M+K'},'mass',{1.007276,22.989218,38.963158}));
      args=processargs(defaults,varargin);
      
      isotopes=Chem.getisotopes(formula,'minabundance',.001);
      if args.dbsave
        msrun=obj.dbsave();
        if isstruct(formula)
          formula=Chem.struct2formula(formula);
        end
        cmd=sprintf('SELECT formula_pk FROM formulas WHERE formula=''%s''',formula);
        fpk=mysql(cmd);
        assert(~isempty(fpk));
      end
      
      tstep=median(diff(obj.time));
      for i=1:length(args.adducts)
        a=args.adducts(i);
        if args.dbsave
          apk=mysql(sprintf('select adduct from adducts where name=''%s''',a.name));
          if isempty(apk)
            fprintf('Added adduct %s to database\n', a.name);
            mysql(sprintf('insert into adducts(name,mass) values(''%s'',%f)',a.name,a.mass));
            apk=mysql(sprintf('select adduct from adducts where name=''%s''',a.name));
            assert(~isempty(apk));
          end
          % Clear any old entries
          mysql(sprintf('DELETE FROM features WHERE msrun=%d AND adduct=%d AND formula_pk=%d',msrun,apk,fpk));
        end
        mz=isotopes(1).mass+a.mass;
        if args.debug
          fprintf('\n%s[%s] m/z=%.4f\n', isotopes(1).name, a.name,isotopes(1).mass+a.mass);
        end
        fl=obj.targetedFeatureDetect(mz,'names',{a.name},'mztol',args.mztol,'noise',args.noise,'timetol',args.timetol);
        fld=fl.deconvolve('noise',500,'oversegmentationfilter',args.timetol);
        for j=1:length(fld.features)
          f=fld.features(j);
          if args.dbsave
            mysql(sprintf('INSERT INTO features(msrun,formula_pk,adduct,rt,mztol) VALUES(%d,%d,%d,%f,%f)',msrun,fpk,apk,f.time,args.mztol));
            feature=mysql('SELECT LAST_INSERT_ID()');
            fprintf('Created feature %d\n', feature);
          end
          if args.debug
            fprintf(' %.2f\n',f.time);
            %fprintf(' %-41.41s %.4f d=%3.0f T=%5.2f, I=%7.0f\n', sprintf('%s@%.2f',isotopes(1).name,f.time), f.mz, (f.mz-mz)*1e4, f.time, f.intensity);
            %fprintf('  timerange=[%.3f,%.3f]\n', f.timerange);
          end
          for k=1:length(isotopes)
            if isotopes(k).abundance*f.intensity < args.noise/2
              continue;
            end
            mzi=isotopes(k).mass+a.mass;
            fli=obj.targetedFeatureDetect(mzi,'names',{isotopes(k).name},'mztol',args.mztol,'noise',0,'timetol',tstep/2,'rt',f.time);
            if args.debug
              fprintf('  %-40.40s %.4f',isotopes(k).name,mzi);
            end
            if length(fli.features)>0
              fi=fli.features(1);
              assert(fi.time==f.time);
              assert(fi.npeaks==1);
              if args.debug
                fprintf(' d=%3.0f',(fi.mz-mzi)*1e4);
              end
              if args.dbsave
                cmd=sprintf(['INSERT INTO isopeaks(feature,isotope,obsmz,ioncount) ',...
                             'SELECT %d,i.isotope,%.5f,%.0f FROM isotopes i ',...
                             'WHERE i.formula=''%s'''],feature, fi.mz,fi.intensity,isotopes(k).name);
                nins=mysql(cmd);
                assert(nins==1);
              end
            else
              if args.debug
                fprintf('      ');
              end
              if args.dbsave
                cmd=sprintf(['INSERT INTO isopeaks(feature,isotope,ioncount) ',...
                             'SELECT %d,i.isotope,0 FROM isotopes i ',...
                             'WHERE i.formula=''%s'''],feature,isotopes(k).name);
                nins=mysql(cmd);
                assert(nins==1);
              end
            end
          end
        end
      end
      if args.debug
        if bestfound>0
          fprintf('Best candidate has %d found: %s mz=%.4f, T=%.2f, \n', bestfound, bestfeature.name, bestmz, bestfeature.time);
        else
          fprintf('No feasible candidates\n');
        end
      end
      if args.debug && bestfound>4
        %keyboard;
      end
    end
      
    function msrun=dbsave(obj)
      if mysql('status')~=0
        mysql('open','35.203.151.202','hts','driver');
      end
      mysql('use hts');
      path=strrep(obj.path,'../','');
      path=strrep(path,'data/','');
      path=strrep(path,'MassSpec/','');
      msrun=mysql(sprintf('select msrun from msruns where filename=''%s''',path));
      if isempty(msrun)
        nins=mysql(sprintf('insert into msruns(created,name,filename) values(NOW(),''%s'',''%s'')',obj.name, path));
        assert(nins==1);
        msrun=mysql(sprintf('select msrun from msruns where filename=''%s''',path));
      end
    end
    
  end % methods
end % classdef