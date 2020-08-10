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
        obj.mzrange=[nanmin(obj.mz),nanmax(obj.mz)];
      end
    end

  end % Methods
end

