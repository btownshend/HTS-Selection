classdef Feature < matlab.mixin.Copyable
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
    snr;    % Only if noise given when constructed
    isotope;   % 1 if this is was identified as the base of a C13 isotope pattern, 2 if next;  0 if not part of a pattern 
  end % properties
  
  methods
    function obj=Feature(peaks,name,noise)
    % Add a feature based on the given peaks (N,[mz,ic,time])
      if nargin<1 || isempty(peaks)
        % Empty constructor
        return;
      end
      [~,maxpeak]=max(peaks(:,2));
      obj.mz=peaks(maxpeak,1);
      obj.time=peaks(maxpeak,3);
      obj.intensity=peaks(maxpeak,2);
      obj.area=sum(peaks(:,2));
      obj.peaks=peaks;
      obj.mzrange=[nanmin(peaks(:,1)),nanmax(peaks(:,1))];
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
      if nargin>=2 && ~isempty(name)
        obj.name=name;
      else
        obj.name=sprintf('%.4f @ %.2f',obj.mz,obj.time);
      end
      if nargin>=3
        obj.snr=obj.intensity/noise;
      else
        obj.snr=nan;
      end
      obj.isotope=0;
    end

    function s=tostring(obj,varargin)
      defaults=struct('mztarget',[],'timetarget',[],'intensitytarget',[]);
      args=processargs(defaults,varargin);
      
    % Return string representation of feature
      s=sprintf('m/z=%.4f',obj.mz);
      if ~isempty(args.mztarget)&& isfinite(args.mztarget)
        dmz=obj.mz-args.mztarget;
        s=sprintf('%s (%3.0f)',s,dmz*10000);
      end
      s=sprintf('%s t=%5.2f',s,obj.time);
      if ~isempty(args.timetarget)&& isfinite(args.timetarget)
        s=sprintf('%s (%5.2f)',s,obj.time-args.timetarget);
      end
      s=sprintf('%s ic=%6.0f',s,obj.intensity);
      if ~isempty(args.intensitytarget) && isfinite(args.intensitytarget)
        s=sprintf('%s (%4.2f)',s,obj.intensity/args.intensitytarget);
      end
      if ~isempty(obj.isotope) && obj.isotope>0
        s=sprintf('%s I%d',s,obj.isotope);
      else
        s=[s,'   '];
      end
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

