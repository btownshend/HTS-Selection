classdef Feature < handle
  properties
    name;
    peaks;   % Peak list
    
    % Computed properties (or imported)
    mz;
    time;
    intensity;
    
    npeaks;
    mzrange;
    timerange;
    fwhh;
    asymmetry;
    tailing;
    area;
  end % properties
  
  methods
    function obj=Feature(peaks,name)
    % Add a feature based on the given peaks (N,[mz,ic,time])
      if nargin<1 || isempty(peaks)
        % Empty constructor
        return;
      end
      if nargin>=2
        obj.name=name;
      else
        obj.name='';
      end
      [~,maxpeak]=max(peaks(:,2));
      obj.mz=peaks(maxpeak,1);
      obj.time=peaks(maxpeak,3);
      obj.intensity=peaks(maxpeak,2);
      obj.area=sum(peaks(:,2));
      obj.peaks=peaks;
      obj.mzrange=[nanmin(peaks(:,1)),nanmax(peaks(:,2))];
      obj.timerange=[min(peaks(:,3)),max(peaks(:,3))];
      p50a=peaks(find(peaks(:,2)>=obj.intensity/2,1),3);
      p50b=peaks(find(peaks(:,2)>=obj.intensity/2,1,'last'),3);
      p10a=peaks(find(peaks(:,2)>=obj.intensity*.1,1),3);
      p10b=peaks(find(peaks(:,2)>=obj.intensity*.1,1,'last'),3);
      p5a=peaks(find(peaks(:,2)>=obj.intensity*.05,1),3);
      p5b=peaks(find(peaks(:,2)>=obj.intensity*.05,1,'last'),3);
      obj.fwhh=p50b-p50a;
      obj.asymmetry=(p10b-obj.time)/(obj.time-p10a);
      obj.tailing=(p5b-p5a)/(2*(obj.time-p5a));
      obj.npeaks=size(peaks,1);
    end

    function set(obj,x)
    % Set from struct x (from csv import)
      fn=fieldnames(x);
      for i=1:length(fn)
        obj.(fn{i})=x.(fn{i});
      end
    end
    
  end % Methods
end

