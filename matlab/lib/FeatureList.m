classdef FeatureList < handle
  properties
    name;
    src;
    params;
    features;   % Array of Feature
  end

  methods(Static)
    function fls=importCSV(filename)
      fls=[];
      x=readtable(filename,'PreserveVariableNames',true);
      % Find peak status row to extract name of this peak list
      peakstatus=find(cellfun(@(z) ~isempty(strfind(z,'Peak status')),x.Properties.VariableNames));
      for i=1:length(peakstatus)
        name=strrep(x.Properties.VariableNames{peakstatus(i)},' Peak status','');
        fl=FeatureList(name,'import');
        s=struct('mz',num2cell(x.([name,' Peak m/z'])),...
                 'name',x.([name,' Peak name']),...
                 'time',num2cell(x.([name,' Peak RT'])),...
                 'intensity',num2cell(x.([name,' Peak height'])),...
                 'npeaks',num2cell(x.([name,' Peak # data points'])),...
                 'mzrange',num2cell([x.([name,' Peak m/z min']),x.([name,' Peak m/z max'])],2),...
                 'timerange',num2cell([x.([name,' Peak RT start']),x.([name,' Peak RT end'])],2),...
                 'fwhh',num2cell(x.([name,' Peak FWHM'])),...
                 'tailing',num2cell(x.([name,' Peak tailing factor'])),...
                 'asymmetry',num2cell(x.([name,' Peak asymmetry factor'])),...
                 'area',num2cell(x.([name,' Peak area'])));
        for j=1:length(s)
          f=Feature();
          f.set(s(j));
          fl.features=[fl.features,f];
        end
        fls=[fls,fl];
      end
    end
    
    function fl=commonfeatures(objs,varargin)
    % Find features that are shared by the lists
      defaults=struct('mztol',.01,'timetol',1,'minmatch',length(objs),'minintensity',1000,'trimtime',7);  
      args=processargs(defaults,varargin);

      args.inputobjs=cellfun(@(z) z.name,objs,'Unif',false);
      fl=FeatureList('','commonfeatures',args);
      all=[];
      for i=1:length(objs)
        f=objs{i}.features;
        sel=[f.intensity]>=args.minintensity & [f.time]>=args.trimtime & [f.time]<=max([f.time])-args.trimtime;
        all=[all,struct('mz',num2cell([f(sel).mz]),'time',num2cell([f(sel).time]),'index',i)];
      end
      fprintf('Have a total of %d features\n', length(all));
      [~,ord]=sort([all.mz]);
      all=all(ord);
      for i=1:length(all)
        if all(i).time<0
          % Already used
          continue;
        end
        for j=i+args.minmatch-1:length(all)
          if all(j).mz-all(i).mz > 2*args.mztol
            break;
          end
          sel=find(abs([all(i+1:j).time]-all(i).time)<args.timetol);
          if length(sel)<args.minmatch-1
            continue;
          end
          aa=all([0,sel]+i);
          if length(unique([aa.index]))>=args.minmatch
            % Found one
            mzrange=[min([aa.mz]),max([aa.mz])];
            timerange=[min([aa.time]),max([aa.time])];
            fprintf('Hit: %d:%d mz=[%.4f,%.4f] t=[%.2f,%.2f]\n', i,j,mzrange,timerange);
            mz=mean(mzrange);
            time=mean(timerange);
            allsel=abs([all.mz]-mz)<args.mztol & abs([all.time]-time)<args.timetol;
            aa=all(allsel);
            if length(unique([aa.index])) < args.minmatch
              keyboard;
            end
            f=Feature();
            f.mz=mz;
            f.time=time;
            % Overload timerange, mzrange with the actual ones by index
            for i=1:length(objs)
              ai=aa([aa.index]==i);
              f.timerange(i)=nanmean([ai.time]);
              f.mzrange(i)=nanmean([ai.mz]);
            end
            fl.features(end+1)=f;
            [all(allsel).time]=deal(-1);  % Blank them out
            break;
          end
        end
      end
      fprintf('Found %d features shared across at least %d/%d input lists\n', length(fl.features), args.minmatch, length(objs));
    end
    
  end
  
  methods
    function obj=FeatureList(name,src,params)
      obj.features=Feature.empty;
      if nargin>0
        obj.name=name;
      end
      if nargin>=2
        obj.src=src;
      end
      if nargin>=3
        obj.params=params;
      end
    end
    
    function r=maptoref(obj,map)
    % Map m/z, time in feature list from 'file' values to 'ref' values 
    % This interpolates using map.mz(:,1) as ref, map.mz(:,2) as file
    % Similar for map.time
      r=FeatureList([obj.name,' toref'],'remap',map);
      r.features=repmat(Feature.empty,size(obj.features));
      for fi=1:length(obj.features)
        f=obj.features(fi);
        peaks=interp1(map.mz(:,2),map.mz(:,1),f.peaks(:,1),'linear','extrap');
        peaks(:,2)=f.peaks(:,2);
        peaks(:,3)=interp1(map.time(:,2),map.time(:,1),f.peaks(:,3),'linear','extrap');
        fnew=Feature(peaks,[f.name,' toref'],f.intensity/f.snr);
        fnew.isotope=f.isotope;
        fnew.labels=f.labels;
        r.features(fi)=fnew;
      end
    end

    function r=maptofile(obj,map)
    % Map m/z, time in feature list from 'ref' values to 'file' values 
    % This interpolates using map.mz(:,1) as ref, map.mz(:,2) as file
    % Similar for map.time
      r=FeatureList([obj.name,' tofile'],'remap',map);
      r.features=repmat(Feature.empty,size(obj.features),f.snr);
      for fi=1:length(obj.features)
        f=obj.features(fi);
        peaks=interp1(map.mz(:,1),map.mz(:,2),f.peaks(:,1),'linear','extrap');
        peaks(:,2)=f.peaks(:,2);
        peaks(:,3)=interp1(map.time(:,1),map.time(:,2),f.peaks(:,3),'linear','extrap');
        fnew=Feature(peaks,[f.name,' tofile'],f.intensity/f.snr);
        fnew.isotope=f.isotope;
        fnew.labels=f.labels;
        r.features(fi)=fnew;
      end
    end

    function append(obj,feature)
      obj.features=[obj.features,feature];
    end
      
    function sortbymz(obj)
      [~,ord]=sort([obj.features.mz]);
      obj.features=obj.features(ord);
    end
    
    function sortbytime(obj)
      [~,ord]=sort([obj.features.time]);
      obj.features=obj.features(ord);
    end
    
    function [fl,sel]=getbymz(obj,mz,varargin)
      defaults=struct('mztol',0.01,'timerange',[]);
      args=processargs(defaults,varargin);
      fl=FeatureList(sprintf('%.4f',mz),'getbymz',args);
      sel=find(abs([obj.features.mz]-mz)<=args.mztol);
      if ~isempty(args.timerange) & ~isempty(sel)
        trange=vertcat(obj.features(sel).timerange);
        % Keep any whose window overlaps the requested window
        sel=sel(~(trange(:,2)>args.timerange(2) | trange(:,1)<args.timerange(1)));
      end
      fl.features=obj.features(sel);
    end
    
    function ploteic(obj,mz,varargin)
    % Add EIC using this feature list to plot, augment legend
      defaults=struct('mztol',0.01,'prefix','','timerange',[-inf,inf]);
      args=processargs(defaults,varargin);

      flmz=obj.getbymz(mz,'mztol',args.mztol,'timerange',args.timerange);
      features=flmz.features;
      if isempty(features)
        fprintf('No features in %s within %.4f of %.4f\n', obj.name, args.mztol, mz);
        return;
      end

      lh=legend();
      leg=get(lh,'String');
      hold on;
      maxdtplot = range([obj.features.time]) * 0.9;  % Avoid plotting wide windows
      for i=1:length(features)
        e=features(i);
        % EIC
        % Only show part of feature that overlaps our time window
        psel=e.peaks(:,3)>=args.timerange(1) & e.peaks(:,3)<=args.timerange(2);
        if any(psel)
          h=plot(e.peaks(psel,3),e.peaks(psel,2));
          leg{end+1}=sprintf('%s.%d %s',args.prefix,find(e==obj.features),e.name);
        else
          h=[];
        end

        % Individual peaks
        ext=e.timerange;
        if diff(ext)<=maxdtplot & ~(ext(1)>args.timerange(2) || ext(2)<args.timerange(1))
          %plot(e.time*[1,1],[0,e.intensity],'Color',get(h(end),'Color'));
          h2=plot([ext(1),e.time,ext(2)],[0.5,1,0.5]*e.intensity,'-','LineWidth',3);
          if isfinite(e.snr)
            plot([ext(1),ext(2)],[1,1]*e.intensity/e.snr,':','HandleVisibility','off');
          end
          if ~isempty(h)
            leg=leg(1:end-1);
            set(h,'HandleVisibility','off');   % No legend entry for above
            set(h2,'Color',get(h,'Color'));
          end
          leg{end+1}=sprintf('%s.%d %s',args.prefix,find(e==obj.features),e.name);
        end
      end
      legend(leg,'Location','best');
    end
    
    function fl=deconvolve(obj,varargin)
      defaults=struct('debug',false,'heightfilter',500,'oversegmentationfilter',0.2,'maxarearatio',33,'maxwidth',2,'trace',0,'minwidth',3.01/60,'noise',500,'minsnr',2);
      args=processargs(defaults,varargin);
      fl=FeatureList([obj.name,' deconvoluted'],'deconvolve',args);
      rejects=struct('accept',0,'snr',0,'area',0,'width',0);
      if args.debug
        fprintf('Deconvolving %d features...',length(obj.features));
      end
      for i=1:length(obj.features)
        if args.debug && mod(i,100)==0
          fprintf('%d...',i);
        end
        p=obj.features(i).peaks;
        %pba=msbackadj(p(:,3),p(:,2),'est','em','smooth','rloess');
        %pbasm=mslowess(p(:,3),pba);
        [pks,pfwhh,pext]=mspeaks(p(:,3),p(:,2),'denoising',true,'heightfilter',args.heightfilter,'oversegmentationfilter',args.oversegmentationfilter,'NoiseEstimator',args.noise,'showplot',ismember(i,args.trace),'style','exttriangle','fwhhfilter',args.minwidth);
        if size(pks,2)==0
          if ismember(i,args.trace)
            keyboard;
          end
          continue;
        end
        minnoise=min(p(p(:,2)>0,2));   % Minimum level acquired
        for j=1:size(pks,1)
          bounds=find(p(:,3)>=pext(j,1) & p(:,3)<=pext(j,2));
          rawpeaks=p(bounds,:);
          rawpeaks=rawpeaks(find(rawpeaks(:,2)>0,1):find(rawpeaks(:,2)>0,1,'last'),:);   % Remove end zeros
          if isempty(rawpeaks)
            % Denoising near a sharp edge occassionally generates a new peak...
            fprintf('mspeaks indicated a peak at %f in feature %d.%d, but no peaks in [%f,%f]\n', pks(j,1), i, j, pext(j,:));
            if ismember(i,args.trace)
              keyboard;
            end
            continue;
          end
          mz=nansum(rawpeaks(:,1).*rawpeaks(:,2))/sum(rawpeaks(:,2));
          % Calculate SNR of peak
          % Use noise estimate algorithm from ADAP paper (not exactly)
          pw=min(100,bounds(end)-bounds(1)+1);
          leftnoise=max(minnoise,p(max(1,bounds(1)-pw*8):bounds(1)-1,2));
          leftnoise=[repmat(minnoise,8*pw-size(leftnoise,1),1);leftnoise];
          rightnoise=max(minnoise,p(bounds(end)+1:min(end,bounds(end)+8*pw),2));
          rightnoise=[rightnoise;repmat(minnoise,8*pw-size(rightnoise,1),1)];
          wind=hanning(2*pw);  wind=wind/sum(wind);
          leftvar =min(conv( leftnoise.^2,wind,'valid'));
          rightvar=min(conv(rightnoise.^2,wind,'valid'));
          noise=sqrt(mean([leftvar;rightvar]));
          f=Feature(rawpeaks,obj.features(i).name,noise);
          if ismember(i,args.trace) || isnan(noise)
            keyboard;
          end
          if f.snr < args.minsnr
            rejects.snr=rejects.snr+1;
            %          elseif f.area/f.intensity > args.maxarearatio
            %rejects.area=rejects.area+1;
          elseif f.fwhh > args.maxwidth
            rejects.width=rejects.width+1;
          else
            fl.append(f);
            %obj.featuresd=[obj.featuresd,struct('mz',mz,'time',pks(j,1),'intensity',pks(j,2),'mzrange',obj.features(i).mzrange,'pfwhh',pfwhh(j,:),'pfext',pext(j,:), 'peaks',rawpeaks,'area',area)];
            rejects.accept=rejects.accept+1;
          end
        end
      end
      if args.debug
        fprintf('%d\n',length(obj.features));
      end
      fl.params.rejects=rejects;
    end

    function al=plotalign(obj,obj2,varargin)
    % Trim time is number of minutes to trim from ends of timewindow
      defaults=struct('mztol',.01,'timetol',1,'minintensity',5000,'trimtime',7);  
      args=processargs(defaults,varargin);

      al=[];
      twindow1=[min([obj.features.time]),max([obj.features.time])];
      f1=obj.features([obj.features.intensity]>=args.minintensity & [obj.features.time]>=twindow1(1)+args.trimtime & [obj.features.time]<=twindow1(2)-args.trimtime  );

      twindow2=[min([obj2.features.time]),max([obj2.features.time])];
      f2=obj2.features([obj2.features.intensity]>=args.minintensity  & [obj2.features.time]>=twindow2(1)+args.trimtime & [obj2.features.time]<=twindow2(2)-args.trimtime  );

      fprintf('Matching %d and %d features...', length(f1), length(f2));
      for i=1:length(f1)
        if mod(i,50)==0
          fprintf('%d...',i);
        end
        f1i=f1(i);
        sel=find(abs([f2.mz]-f1i.mz)<=args.mztol & abs([f2.time]-f1i.time)<=args.timetol);
        if length(sel)>0
          al=[al;[repmat(f1i.mz,length(sel),1),[f2(sel).mz]',repmat(f1i.time,length(sel),1),[f2(sel).time]']];
        end
      end
      fprintf('\n');
      ti=sprintf('%s vs %s',obj.name,obj2.name);
      setfig(ti);clf;
      layout=tiledlayout('flow');
      nexttile
      plot(al(:,1),al(:,2)-al(:,1),'o');
      xlabel([obj.name,' m/z']);
      ylabel('\Delta m/z');
      nexttile
      tquant=(1/60)/2;
      plot(al(:,3)+(rand(size(al(:,3)))-0.5)*tquant,al(:,4)-al(:,3)+(rand(size(al(:,3)))-0.5)*tquant,'o');
      xlabel([obj.name,' time']);
      ylabel('\Delta time');
      title(layout,ti);
    end
    
    function ilist=findisotopes(obj,varargin)
      defaults=struct('mztol',0.004,'timetol',0.1,'minratio',1e-4,'maxratio',13e-4,'debug',false);
      args=processargs(defaults,varargin);

      ilist=[];
      isodiff=1.0030;   % Although C13-C12 should be 1.003355, we see 1.0030 empircally
      % cratio=.0107/.9893;
      % ncarbons=[8,30];
      % ratio1=binopdf(1,ncarbons,cratio);
      % ratio2=binopdf(2,ncarbons,cratio);
      % fprintf('For [%s] carbons expect ratios of [%s], [%s]\n', sprintf('%d ',ncarbons), sprintf('%.3f ',ratio1), sprintf('%.3f ',ratio2./ratio1));
      allmz=[obj.features.mz];
      for i=1:length(obj.features)
        fi=obj.features(i);
        mzsel=find(abs(allmz-fi.mz-isodiff) < args.mztol);
        for jj=1:length(mzsel)
          j=mzsel(jj);
          fj=obj.features(j);
          if abs(fj.time-fi.time) < args.timetol
            iratio=fj.intensity/fi.intensity;
            if iratio>=args.minratio*fi.mz && iratio<=args.maxratio*fi.mz
              ilist(end+1,:)=[i,j,fi.mz,fj.mz-fi.mz,fj.time-fi.time,iratio];
            end
          end
        end
      end

      % Add a column to indicate the ordering of a sequence of isotopes
      order=zeros(length(obj.features),1);
      order(setdiff(ilist(:,1),ilist(:,2)))=1;   % If it appears as a base, but never as a secondary, then assign it order 1
      for i=1:10
        sel=ismember(ilist(:,1),find(order==i));
        if ~any(sel)
          break;
        end
        order(ilist(sel,2))=i+1;   % assign to secondaries
      end

      % Set in features
      tmp=num2cell(order);
      [obj.features.isotope]=tmp{:};
      %for i=1:length(obj.features)
      %        obj.features(i).isotope=order(i);
      % end
      ilist(:,end+1)=order(ilist(:,1));
      ilist(:,end+1)=order(ilist(:,2));

      fprintf('%d/%d (%.0f%%) features appear to be in isotope patterns\n',sum(order>0),length(order),mean(order>0)*100);
    end

    function loadlabels(obj, filename, varargin)
    % Load and apply labels from CSV file
      defaults=struct('mztol',0.004,'timetol',0.1,'debug',false,'minintensity',5000);
      args=processargs(defaults,varargin);

      x=readtable(filename);
      if ~ismember('exact_adduct_mass',x.Properties.VariableNames)
        x=xrenamevars(x,'mz','exact_adduct_mass');
      end
      if ~all(ismember({'id','name','exact_adduct_mass','ion_form'},x.Properties.VariableNames))
        error('Missing some required variables in %s\n', filename);
      end
      if ~ismember('RT',x.Properties.VariableNames)
        x.RT(:)=0;
      end
      for i=1:length(obj.features)
        f=obj.features(i);
        if f.intensity<args.minintensity
          continue;
        end
        sel=find(abs(f.mz-x.exact_adduct_mass)<args.mztol & (x.RT==0 | abs(f.time-x.RT)<args.timetol));
        [a,b,c]=unique(x.id(sel));
        sel=sel(b);
        if ~isempty(sel)
          fprintf('Matched feature %d %.4f@%.2f to:',i,f.mz,f.time);
          for jj=1:length(sel)
            j=sel(jj);
            fprintf(' %s%s %.4f;',x.name{j}, x.ion_form{j}, x.exact_adduct_mass(j));
            if x.RT(j)~=0
              fprintf('@%.2f',x.RT(j));
            end
          end
          fprintf('\n');
          if length(unique(x.name(sel))) ~= length(sel)
            keyboard;
          end
        end
      end
    end
    
  end % methods
end % classdef