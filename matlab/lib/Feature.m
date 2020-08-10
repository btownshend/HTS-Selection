classdef Feature < handle
  properties
    mz;
    time;
    intensity;
    mzrange;
    area;
    peaks;
  end % properties
  
  methods
    function obj=Feature(peaks,mzrange)
    % Add a feature based on the given peaks (N,[mz,ic,time])
      [~,maxpeak]=max(peaks(:,2));
      obj.mz=peaks(maxpeak,1);
      obj.time=peaks(maxpeak,3);
      obj.intensity=peaks(maxpeak,2);
      obj.area=sum(peaks(:,2));
      obj.peaks=peaks;
      if nargin>=2 && ~isempty(mzrange)
        obj.mzrange=mzrange;
      else
        obj.mzrange=[nanmin(peaks(:,1)),nanmax(peaks(:,2))];
      end
    end

    function t=gettimerange(obj)
      t=[min(obj.peaks(:,3)),max(obj.peaks(:,3))];
    end
    
    function t=getfwhh(obj)
      t=obj.peaks([find(obj.peaks(:,2)>=obj.intensity/2,1),find(obj.peaks(:,2)>=obj.intensity/2,1,'last')],3)';
    end
    
  end % Methods
end

