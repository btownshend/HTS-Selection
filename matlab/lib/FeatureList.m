classdef FeatureList < handle
  properties
    name;
    features;   % Array of Feature
  end
  
  methods
    function obj=FeatureList(name)
      obj.name=name;
      obj.features=[];
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
      defaults=struct('debug',false,'heightfilter',5000,'oversegmentationfilter',60,'maxarearatio',33,'maxwidth',120,'trace',0);
      % Note: when oversegmentationfilter is active, the fwhh includes the merged peak's fwhh even if they are much lower than the main one
      args=processargs(defaults,varargin);
      fl=FeatureList([obj.name,' deconvoluted']);
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
          area=sum(rawpeaks(:,2));
          mz=nansum(rawpeaks(:,1).*rawpeaks(:,2))/sum(rawpeaks(:,2));
          if ismember(i,args.trace)
            keyboard;
          end
          if area/pks(j,2) < args.maxarearatio & diff(pfwhh(j,:))<=args.maxwidth
            f=Feature(rawpeaks);
            fl.append(f);
            %obj.featuresd=[obj.featuresd,struct('mz',mz,'time',pks(j,1),'intensity',pks(j,2),'mzrange',obj.features(i).mzrange,'pfwhh',pfwhh(j,:),'pfext',pext(j,:), 'peaks',rawpeaks,'area',area)];
          end
        end
      end
    end
    
    function importCSV(obj,filename)
      x=readtable(filename);
    end
    
  end % methods
end % classdef