classdef Feature < handle
  properties
    mz;
    time;
    intensity;
    mzrange;
    timerange;
    fwhh;
    asymmetry;
    tailing;
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
      obj.timerange=[min(peaks(:,3)),max(peaks(:,3))];
      p50a=peaks(find(peaks(:,2)>=obj.intensity/2,1),3);
      p50b=peaks(find(peaks(:,2)>=obj.intensity/2,1,'last'),3);
      p10a=peaks(find(peaks(:,2)>=obj.intensity*.1,1),3);
      p10b=peaks(find(peaks(:,2)>=obj.intensity*.1,1,'last'),3);
      p5a=peaks(find(peaks(:,2)>=obj.intensity*.05,1),3);
      p5b=peaks(find(peaks(:,2)>=obj.intensity*.05,1,'last'),3);
      obj.fwhh=p50b-p50a;
      obj.asymmetry=(obj.time-p10a)/(p10b-obj.time);
      obj.tailing=(obj.time-p5a)*(p5b-p5a)/(2*(obj.time-p5a)*(p5b-obj.time));
      obj.npeaks=size(peaks,1);
    end

    end
    
  end % Methods
end

