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
      for fi=1:length(obj.features)
        f=obj.features(fi);
        peaks=interp1(map.mz(:,2),map.mz(:,1),f.peaks(:,1),'linear','extrap');
        peaks(:,2)=f.peaks(:,2);
        peaks(:,3)=interp1(map.time(:,2),map.time(:,1),f.peaks(:,3),'linear','extrap');
        r.append(Feature(peaks));
      end
    end

    function r=maptofile(obj,map)
    % Map m/z, time in feature list from 'ref' values to 'file' values 
    % This interpolates using map.mz(:,1) as ref, map.mz(:,2) as file
    % Similar for map.time
      r=FeatureList([obj.name,' tofile'],'remap',map);
      for fi=1:length(obj.features)
        f=obj.features(fi);
        peaks=interp1(map.mz(:,1),map.mz(:,2),f.peaks(:,1),'linear','extrap');
        peaks(:,2)=f.peaks(:,2);
        peaks(:,3)=interp1(map.time(:,1),map.time(:,2),f.peaks(:,3),'linear','extrap');
        r.append(Feature(peaks));
      end
    end

    function append(obj,feature)
      obj.features=[obj.features,feature];
    end
      
    function sortbymz(obj)
      [~,ord]=sort([obj.features.mz]);
      obj.features=obj.features(ord);
    end
    
    function [fl,sel]=getbymz(obj,mz,varargin)
      defaults=struct('mztol',0.01);
      args=processargs(defaults,varargin);
      fl=FeatureList(sprintf('%.4f',mz),'getbymz',args);
      sel=find(abs([obj.features.mz]-mz)<=args.mztol);
      fl.features=obj.features(sel);
    end
    
    function ploteic(obj,mz,varargin)
    % Add EIC using this feature list to plot, augment legend
      defaults=struct('mztol',0.01,'prefix','','timerange',[-inf,inf]);
      args=processargs(defaults,varargin);

      flmz=obj.getbymz(mz,'mztol',args.mztol);
      features=flmz.features;
      if isempty(features)
        fprintf('No features in %s within %.4f of %.4f\n', obj.name, args.mztol, mz);
        return;
      end

      lh=legend();
      leg=get(lh,'String');
      hold on;
      for i=1:length(features)
        e=features(i);
        % EIC
        psel=e.peaks(:,3)>=args.timerange(1) & e.peaks(:,3)<=args.timerange(2);
        if any(psel)
          h=plot(e.peaks(psel,3),e.peaks(psel,2));
        else
          h=[];
        end

        % Individual peaks
        ext=e.timerange;
        if diff(ext)<30 & (ext(1)<args.timerange(2) || ext(2)>args.timerange(1))
          %plot(e.time*[1,1],[0,e.intensity],'Color',get(h(end),'Color'));
          h2=plot([ext(1),e.time,ext(2)],[0.5,1,0.5]*e.intensity,'-','LineWidth',3);
          if ~isempty(h)
            set(h,'HandleVisibility','off');   % No legend entry for above
            set(h2,'Color',get(h,'Color'));
          end
        end
        leg{end+1}=sprintf('%s.%d %s',args.prefix,find(e==obj.features),e.name);
      end
      legend(leg,'Location','best');
    end
    
    function fl=deconvolve(obj,varargin)
      defaults=struct('debug',false,'heightfilter',1000,'oversegmentationfilter',0.3,'maxarearatio',33,'maxwidth',120/60,'trace',0,'minwidth',0.2);
      args=processargs(defaults,varargin);
      fl=FeatureList([obj.name,' deconvoluted'],'deconvolve',args);
      for i=1:length(obj.features)
        %fprintf('%d ',i);
        p=obj.features(i).peaks;
        %pba=msbackadj(p(:,3),p(:,2),'est','em','smooth','rloess');
        %pbasm=mslowess(p(:,3),pba);
        [pks,pfwhh,pext]=mspeaks(p(:,3),p(:,2),'denoising',true,'heightfilter',args.heightfilter,'oversegmentationfilter',args.oversegmentationfilter,'showplot',ismember(i,args.trace),'style','fwhhtriangle','fwhhfilter',args.minwidth);
        if size(pks,2)==0
          continue;
        end
        for j=1:size(pks,1)
          rawpeaks=p(p(:,3)>=pext(j,1) & p(:,3)<=pext(j,2),:);
          rawpeaks=rawpeaks(find(rawpeaks(:,2)>0,1):find(rawpeaks(:,2)>0,1,'last'),:);   % Remove end zeros
          mz=nansum(rawpeaks(:,1).*rawpeaks(:,2))/sum(rawpeaks(:,2));
          if ismember(i,args.trace)
            keyboard;
          end
          f=Feature(rawpeaks);
          if f.area/f.intensity < args.maxarearatio & f.fwhh<=args.maxwidth
            fl.append(f);
            %obj.featuresd=[obj.featuresd,struct('mz',mz,'time',pks(j,1),'intensity',pks(j,2),'mzrange',obj.features(i).mzrange,'pfwhh',pfwhh(j,:),'pfext',pext(j,:), 'peaks',rawpeaks,'area',area)];
          end
        end
      end
    end
    
  end % methods
end % classdef