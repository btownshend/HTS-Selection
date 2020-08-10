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
    
  end % methods
end % classdef