classdef MassSpec < handle
  properties
    path;    % Full path
    name;    % Short name for this file
    moles;   % Nominal number of moles (per compound) loaded
    peaks;  % Peaks from mzxml2peaks;  peaks{i}(k,1:2) is [mz,ic] for elution i, at peak k 
    time;   % Time of elutions (in minutes)
    mzrange;   % Range [low,high] of m/z to plot/consider
    resamp;   % struct (mz,y,n) of uniformly sampled data
    clusters;	% Cluster peaks (across m/z and time)
    clustersettings;
    ident;    % struct array of identified peaks
    featurelists;   % Array of feature lists
    mzxml;
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
      for i=1:length(c)
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
      for i=1:args.nlist
        fprintf('%8.4f %s\n',c(i,1),sprintf('%4.2f ',c(i,2:end)/c(i,args.ref+1)));
      end
    end
  end
  
  methods
    function obj=MassSpec(path,varargin)
      defaults=struct('debug',false);
      args=processargs(defaults,varargin);
      fprintf('Loading %s...\n', path);
      mzxml=mzxmlread(path);
      [obj.peaks,obj.time]=mzxml2peaks(mzxml);
      obj.time=obj.time/60;   % Convert to minutes
      obj.mzrange=[min(cellfun(@(z) z(1,1), obj.peaks)),max(cellfun(@(z) z(end,1), obj.peaks))];
      obj.path=path;
      z=strsplit(obj.path,'/');
      obj.name=z{end};
      obj.mzxml=mzxml;
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

    function clusterpeaks(obj,varargin)
    % Cluster all significant peaks and return a table of results
    % Used for de novo analysis (without info on expected components)
      defaults=struct('mzwindow',0.01,'timewindow',200/60, 'minic',1e3, 'noise',300,'maxpeaks',2000,'noisefrac',.01);
      args=processargs(defaults,varargin);

      if args.noise>args.minic
        args.noise=args.minic;
      end
      
      peaks=obj.peaks;   % Copy peaks so we can clear it
      res=[];
      for pnum=1:args.maxpeaks
        maxp=zeros(length(peaks),2);
        for i=1:length(peaks)
          peak=peaks{i};
          [maxp(i,1),maxp(i,2)]=max(peak(:,2));
        end
        [maxmax,maxxpos]=max(maxp(:,1));
        if isempty(maxmax) || maxmax<args.minic
          break;
        end
        mz=peaks{maxxpos}(maxp(maxxpos,2),1);
        trange=find(abs(obj.time-obj.time(maxxpos))<=args.timewindow);
        ic=zeros(size(trange));
        for i=1:length(trange)
          mzrange=find(abs(peaks{trange(i)}(:,1)-mz)<args.mzwindow);
          if ~isempty(mzrange)
            ic(i)=sum(peaks{trange(i)}(mzrange,2));
            peaks{trange(i)}(mzrange,2)=0;
          end
        end
        thresh=max([args.noise,args.noisefrac*maxmax]);
        sel=min(find(ic>thresh)):max(find(ic>thresh));
        tvalid=trange(sel);
        ic=ic(sel);
        ictotal=sum(ic);
        res=[res,struct('mz',mz,'elution',maxxpos,'elutionrange',[min(tvalid),max(tvalid)],'maxic',maxmax,'ic',ic,'ictotal',ictotal,'thresh',thresh)];
        % fprintf('Max peak of %7.0f (Q=%.2f) at m/z=%.4f, t=%.2f [%.2f,%.2f]\n', ictotal, maxmax/ictotal, mz, obj.time([maxxpos,min(tvalid),max(tvalid)]));
      end
      [~,ord]=sort([res.ictotal],'desc');
      res=res(ord);
      obj.clusters=res;
      obj.clustersettings=args;
      fprintf('Found a total of %d/%d peaks with IC>=%.0f>=%.0f\n', length(res),args.maxpeaks, min([res.maxic]), args.minic);
    end
    
    function [ic,mz,time]=mzscan(obj, mztarget, varargin)
    % Get ion counts for given mztarget across elution times
    % Also return weighted mean m/z for each elution time
      defaults=struct('mztol',0.01,'timerange',[]);
      args=processargs(defaults,varargin);
      % Build a table of elution vs IC for this M/Z
      if isempty(args.timerange)
        first=1;
        last=length(obj.peaks);
      else
        first=find(obj.time>=args.timerange(1),1);
        last=find(obj.time<=args.timerange(2),1,'last');
      end
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
    % Find compound directly (not using clusters)
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

    function listidents(obj)
    %  List all the identified compounds
      for i=1:length(obj.ident)
        id=obj.ident(i);
        fprintf('%2d %-12.12s M/Z=%8.4f\n',i, id.name,id.mztarget);
        for k=1:length(id.peaks)
          pk=id.peaks(k);
          fprintf('   %12.12s     %8.4f  T=[%6.2f,%6.2f,%6.2f]  IC=%8.0f\n',  '',pk.mzobs, pk.elutions, pk.ictotal);
        end
      end
    end

    function plotident(obj,i)
    % Make plot for i-th identification entry
      id=obj.ident(i);
      if isempty(id.peaks)
        fprintf('No peaks for ident(%d)\n',i);
        return;
      end
      
      setfig(sprintf(' %.4f',id.mztarget));clf;

      subplot(211);
      pk=id.peaks(1);
      semilogy(obj.time,id.ioncount);
      hold on;
      ax=axis;
      tmargin=ceil((pk.elutions(3)-pk.elutions(1))/4);
      ax(1)=obj.time(max(1,pk.elutions(1)-tmargin));
      ax(2)=obj.time(min(length(obj.peaks),pk.elutions(3)+tmargin));
      axis(ax);
      plot(ax(1:2),pk.icthresh*[1,1],':r');
      plot(obj.time(pk.elutions(1))*[1,1],ax(3:4),':r');
      plot(obj.time(pk.elutions(2))*[1,1],ax(3:4),':g');
      plot(obj.time(pk.elutions(3))*[1,1],ax(3:4),':r');
      xlabel('Time (min)');
      ylabel('Ion count');
      title(sprintf('M/Z %.4f (err=%.4f) IC=%.0f',id.mztarget,pk.mzobs-id.mztarget,pk.ictotal));

      subplot(212);
      stem(obj.peaks{pk.elutions(2)}(:,1),obj.peaks{pk.elutions(2)}(:,2));
      xlabel('m/z');
      ylabel('Ion count');
      title(sprintf('T(%d)=%f',pk.elutions(2),obj.time(pk.elutions(2))));
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
      defaults=struct('mztol',0.01,'debug',false,'mingroupsize',5,'groupthresh',500,'minintensity',1000,'noise',500);
      args=processargs(defaults,varargin);
      % Build list of all peaks with time index as 3rd column
      allpks=[];
      for i=1:length(obj.peaks)
        p=obj.peaks{i};
        p=p(p(:,2)>=args.noise,:);
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
      fl=FeatureList([obj.name,' chromatograms']);
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
            fprintf('m/z=%.4f, IC=%.0f new EIC [%.4f,%.4f] with %d scans (%d>%.0f) pass=%d\n', ...
                    allpks(i,1:2), mzrange, length(sel),length(gscans),args.groupthresh, ok);
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

    function deconvolve(obj,varargin)
      fl=obj.featurelists(end).deconvolve(varargin{:});
      obj.featurelists=[obj.featurelists,fl];
    end      

    function ploteic(obj,mz,varargin)
      defaults=struct('debug',false,'mztol',0.01);
      args=processargs(defaults,varargin);

      if isempty(obj.featurelists)
        error('No features -- need to buildchromatograms()');
      elseif length(obj.featurelists)<2
        fprintf('Warning - no deconvolved chromatograms\n');
      end
      
      ti=sprintf('Chromatogram %.4f - %.4f',mz-args.mztol,mz+args.mztol);
      setfig(ti);clf;
      leg={};
      h=[];
      % Plot base traces from raw peaks
      [ic,~,time]=obj.mzscan(mz,'mztol',args.mztol);
      h(end+1)=plot(time,ic);
      leg{end+1}='Base peaks';
      hold on;
      
      %sel=find(arrayfun(@(z) z.mzrange(1)-args.mztol<=mz && z.mzrange(2)+args.mztol>=mz,obj.eic));
      for ifl=1:length(obj.featurelists)
        fl=obj.featurelists(ifl);
        features=fl.getbymz(mz,'mztol',args.mztol);
        if isempty(features)
          fprintf('No features in %s within %.4f of %.4f\n', fl.name, args.mztol, mz);
        end

        for i=1:length(features)
          e=features(i);
          if size(e.peaks,1) == length(obj.time)
            % EIC
            h(end+1)=plot(e.peaks(:,3),e.peaks(:,2));
            leg{end+1}=e.name;
          else
            % Deconvolved EIC
            ext=e.timerange;
            plot(e.time*[1,1],[0,e.intensity],'Color',get(h(end),'Color'));
            h(end+1)=plot([ext(1),e.time,ext(2)],[0.5,1,0.5]*e.intensity,'-','LineWidth',3);
            leg{end+1}=e.name;
          end
        end
      end
      xlabel('Time');
      ylabel('Intensity');
      legend(h,leg);
      title(ti);
    end
    
    function plotmz(obj,time)
    % Plot m/z vs ioncount for trace nearest given time
      if isempty(obj.resamp)
        obj.resample();
      end
      [~,closest]=min(abs(time-obj.time));
      ti=sprintf('m/z @ T=%.2f',obj.time(closest));
      setfig(ti);clf;
      plot(obj.resamp.mz,obj.resamp.y(:,closest));
      title(ti)
      xlabel('m/z')
      ylabel('Relative Intensity')
      hold on;
      ax=axis;
      for i=1:length(obj.ident)
        id=obj.ident(i);
        plot(id.mztarget*[1,1],[0,ax(4)*1.1],':r')
      end
      % List major peaks
      pks=obj.peaks{closest};
      [~,ord]=sort(pks(:,2),'desc');
      for ii=1:10
        i=ord(ii);
        fprintf('M/Z=%.4f, IC=%6.0f\n', pks(i,:));
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
      for i=1:length(obj.ident)
        id=obj.ident(i);
        for k=1:length(id.peaks)
          pk=id.peaks(k);
          if k==1
            col='g';  % Green for main peak
          else
            col='r';   % Red for others
          end
          plot(obj.time(pk.elutions(2))*[1,1],[0,ax(4)*1.1],[':',col])
        end
      end
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
    
    
    function x=compare(obj,o2)
    % Compare identified peaks with another MassSpec object
      [c,ia,ib]=intersect({obj.ident.name},{o2.ident.name});

      fprintf('Matched %d of %d+%d identified compounds\n', length(c), length(obj.ident), length(o2.ident));
      fprintf('                         %20.20s                      %20.20s\n', obj.name, o2.name);
      %fprintf('Name           M/Z Exp.   dM/Z Elution Ion Cnt        dM/Z Elution Ion Cnt  IC 2/1\n');
      elutions=[];
      x=table();
      for i=1:length(c)
        i1=obj.ident(ia(i));
        i2=o2.ident(ib(i));
        if i1.mztarget~=i2.mztarget
          fprintf('Warning: %s has different M/Z targets (%.3f, %.3f)\n', i1.mztarget, i2.mztarget);
        end

        peakdelta=1e10;
        pk1=1;pk2=1;
        for k1=1:length(i1.peaks)
          t1=obj.time(i1.peaks(k1).elutions(2));
          for k2=1:length(i2.peaks)
            t2=o2.time(i2.peaks(k2).elutions(2));
            if abs(t1-t2)<abs(peakdelta)
              pk1=k1;
              pk2=k2;
              peakdelta=t2-t1;
            end
          end
        end
        if peakdelta==1e10
          continue;
        end
        p1=i1.peaks(pk1);
        p2=i2.peaks(pk2);
        row=table();
        row.name=i1.name;
        row.mztarget=i1.mztarget;
        row.id1=sprintf('%4d.%-3d',ia(i),pk1);
        row.mzdelta1=round(p1.mzobs-i1.mztarget,3);
        row.elution1=obj.time(p1.elutions(2));
        row.ic1=round(p1.ictotal,0);

        row.id2=sprintf('%4d.%-2d',ib(i),pk2);
        row.mzdelta2=round(p2.mzobs-i1.mztarget,3);
        row.elution2=o2.time(p2.elutions(2));
        row.ic2=round(p2.ictotal,0);
        row.tdelta=round(peakdelta,2);
        row.flags='  ';
        if abs(p1.mzobs-p2.mzobs)>0.005
          row.flags(1)='M';
        end
        if abs(peakdelta)>10
          row.flags(2)='T';
        end
        %fprintf('%-12.12s   %7.3f  %6.3f %7.2f %7.0f      %6.3f %7.2f %7.0f %6.1f%%\n', c{i}, i1.mztarget, p1.mzobs-i1.mztarget, p1.elutions(2), p1.ictotal, i2.mzobs-i1.mztarget, p2.elutions(2), p2.ictotal, p2.ictotal/p1.ictotal*scaling*100);
        x=vertcat(x,row);
      end
      scaling=median(x.ic1./x.ic2);   % Scale factor to match ion counts
      x.icratio=round(100*x.ic2./x.ic1*scaling,1);
      disp(x)
      ti=sprintf('%s vs %s',obj.name, o2.name);

      setfig(ti);clf;
      subplot(221);
      loglog(x.ic1, x.ic2,'o');
      xlabel(sprintf('Ion Count %s',obj.name));
      ylabel(sprintf('Ion Count %s',o2.name));
      ax=axis;
      hold on;
      plot(ax(1:2),ax(1:2)/scaling,':r');
      title(ti);
      subplot(222);
      plot(x.mztarget, x.mzdelta1,'o');
      hold on;
      plot(x.mztarget, x.mzdelta2,'x');
      xlabel('M/Z Exact')
      ylabel('M/Z Diff')
      subplot(223);
      plot(x.elution1,x.elution2,'o');
      xlabel('T1 (sec)');
      ylabel('T2 (sec)')
    end
    
    function ilist=findisotopes(obj,pknum,varargin)
      defaults=struct('mztol',0.01,'debug',false);
      args=processargs(defaults,varargin);

      if nargin<2 || isempty(pknum)
        ilist=[];
        for i=1:length(obj.peaks)
          ilist=[ilist;obj.findisotopes(i,'mztol',args.mztol,'debug',args.debug)];
        end
        return;
      end
      isodiff=1.003355;   % C13 - C12 mass
      pks=obj.peaks{pknum};
      ilist=[];
      i1=find(pks(:,1)>obj.mzrange(1),1);
      i2=find(pks(:,1)<=obj.mzrange(2),1,'last');
      for i=i1:i2
        for j=i+1:i2
          if abs(pks(j,1)-pks(i,1)-isodiff)<args.mztol
            ilist(end+1,:)=[pknum,pks(i,:),pks(j,:)];
            if args.debug
              fprintf('Trace %d, Isotope: %.4f +%.4f  %.1f%%\n', pknum, pks(i,1), pks(j,1)-pks(i,1), pks(j,2)/pks(i,2)*100);
            end
          elseif pks(j,1)-pks(i,1)-isodiff>args.mztol
            break;
          end
        end
      end
    end
    
  end % methods
end % classdef