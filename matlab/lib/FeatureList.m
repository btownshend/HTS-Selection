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
      obj.name=name;
      obj.features=[];
      if nargin>=2
        obj.src=src;
      end
      if nargin>=3
        obj.params=params;
      end
    end
    
    function append(obj,feature)
      obj.features=[obj.features,feature];
    end
      
    function sortbymz(obj)
      [~,ord]=sort([obj.features.mz]);
      obj.features=obj.features(ord);
    end
    
    function f=getbymz(obj,mz,varargin)
      defaults=struct('mztol',0.01);
      args=processargs(defaults,varargin);
      f=obj.features(abs([obj.features.mz]-mz)<=args.mztol);
    end
    
    function fl=deconvolve(obj,varargin)
      defaults=struct('debug',false,'heightfilter',1000,'oversegmentationfilter',60/60,'maxarearatio',33,'maxwidth',120/60,'trace',0);
      % Note: when oversegmentationfilter is active, the fwhh includes the merged peak's fwhh even if they are much lower than the main one
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