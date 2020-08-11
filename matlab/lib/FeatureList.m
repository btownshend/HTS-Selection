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
    
    function fl=deconvolve(obj,varargin)
      defaults=struct('debug',false,'heightfilter',1000,'oversegmentationfilter',60/60,'maxarearatio',33,'maxwidth',120/60,'trace',0);
      args=processargs(defaults,varargin);
      fl=FeatureList([obj.name,' deconvoluted'],'deconvolve',args);
      for i=1:length(obj.features)
        %fprintf('%d ',i);
        p=obj.features(i).peaks;
        %pba=msbackadj(p(:,3),p(:,2),'est','em','smooth','rloess');
        %pbasm=mslowess(p(:,3),pba);
        [pks,pfwhh,pext]=mspeaks(p(:,3),p(:,2),'denoising',true,'heightfilter',args.heightfilter,'oversegmentationfilter',args.oversegmentationfilter,'showplot',ismember(i,args.trace),'style','fwhhtriangle');
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