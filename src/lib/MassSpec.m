classdef MassSpec < handle
  properties
    path;    % Full path
    name;    % Short name for this file
    moles;   % Nominal number of moles (per compound) loaded
    peaks;  % Peaks from mzxml2peaks;  peaks{i}(k,1:2) is [mz,time] for elution i, at peak k 
    time;   % Time of elutions
    mzrange;   % Range [low,high] of m/z to plot/consider
    resamp;   % struct (mz,y,n) of uniformly sampled data
    clusters;	% Cluster peaks (across m/z and time)
    clustersettings;
    ident;    % struct array of identified peaks
    mzoffset;   % M/Z offset (obs-actual) -- for info only, already applied to the other fields
  end
  
  properties(Constant)
    mzres=0.01;
  end
  
  methods
    function obj=MassSpec(path,varargin)
      defaults=struct('mzoffset',0.0);  % mzoffset is observed-actual M/Z
      args=processargs(defaults,varargin);
      fprintf('Loading %s...\n', path);
      mzxml=mzxmlread(path);
      [obj.peaks,obj.time]=mzxml2peaks(mzxml);
      obj.mzrange=[min(cellfun(@(z) z(1,1), obj.peaks)),max(cellfun(@(z) z(end,1), obj.peaks))];
      obj.mzoffset=0;
      if args.mzoffset~=0
        obj.adjmzoffset(args.mzoffset);
      end
      obj.path=path;
      z=strsplit(obj.path,'/');
      obj.name=z{end};
    end

    function adjmzoffset(obj,mzoffset)
    % Correct the M/Z
      for i=1:length(obj.peaks)
        obj.peaks{i}(:,1)=obj.peaks{i}(:,1)-(mzoffset-obj.mzoffset);   % Adjust offset
      end
      obj.mzoffset=mzoffset;
    end
      
    function setLoad(obj,moles)
      obj.moles=moles;
    end
    
    function filter(obj,trange,mzrange)
      if nargin>=2 && ~isempty(trange)
        sel=find(obj.time>=trange(1) & obj.time<=trange(2));
        obj.time=obj.time(sel);
        obj.peaks=obj.peaks(sel);
      end
      if nargin>=3 && ~isempty(mzrange)
        obj.mzrange=mzrange;
      end
    end

    function resample(obj)
      n=round(diff(obj.mzrange)/obj.mzres+1);
      fprintf('Resampling peaks...');
      [mz,y]=msppresample(obj.peaks,n,'Range',obj.mzrange,'FWHH',obj.mzres);
      fprintf('done\n');
      obj.resamp=struct('mz',mz,'y',y,'n',n);
    end

    function heatmap(obj,markers)
      if isempty(obj.resamp)
        obj.resample();
      end
      if nargin<2
        markers=[];
      end
      if isempty(markers)
        msheatmap(obj.resamp.mz,obj.time,log(obj.resamp.y),'resolution',obj.mzres);
      else
        msheatmap(obj.resamp.mz,obj.time,log(obj.resamp.y),'resolution',obj.mzres,'markers',markers);
      end
    end

    function clusterpeaks(obj,varargin)
    % Cluster all significant peaks and return a table of results
    % Used for de novo analysis (without info on expected components)
      defaults=struct('mzwindow',0.01,'timewindow',200, 'minic',1e3, 'noise',300,'maxpeaks',2000,'noisefrac',.01);
      args=processargs(defaults,varargin);

      if args.noise>args.minic
        args.noise=args.minic;
      end
      
      peaks=obj.peaks;   % Copy peaks so we can clear it
      res=[];
      for pnum=1:args.maxpeaks
        maxp=zeros(length(peaks),2);
        for i=1:length(peaks)
          peak=peaks{i};
          [maxp(i,1),maxp(i,2)]=max(peak(:,2));
        end
        [maxmax,maxxpos]=max(maxp(:,1));
        if isempty(maxmax) || maxmax<args.minic
          break;
        end
        mz=peaks{maxxpos}(maxp(maxxpos,2),1);
        trange=find(abs(obj.time-obj.time(maxxpos))<=args.timewindow);
        ic=zeros(size(trange));
        for i=1:length(trange)
          mzrange=find(abs(peaks{trange(i)}(:,1)-mz)<args.mzwindow);
          if ~isempty(mzrange)
            ic(i)=sum(peaks{trange(i)}(mzrange,2));
            peaks{trange(i)}(mzrange,2)=0;
          end
        end
        thresh=max([args.noise,args.noisefrac*maxmax]);
        sel=min(find(ic>thresh)):max(find(ic>thresh));
        tvalid=trange(sel);
        ic=ic(sel);
        ictotal=sum(ic);
        res=[res,struct('mz',mz,'elution',maxxpos,'elutionrange',[min(tvalid),max(tvalid)],'maxic',maxmax,'ic',ic,'ictotal',ictotal,'thresh',thresh)];
        % fprintf('Max peak of %7.0f (Q=%.2f) at m/z=%.4f, t=%4.0f [%4.0f,%4.0f]\n', ictotal, maxmax/ictotal, mz, obj.time([maxxpos,min(tvalid),max(tvalid)]));
      end
      [~,ord]=sort([res.ictotal],'desc');
      res=res(ord);
      obj.clusters=res;
      obj.clustersettings=args;
      fprintf('Found a total of %d/%d peaks with IC>=%.0f>=%.0f\n', length(res),args.maxpeaks, args.minic,min([res.maxic]));
    end
    
    function [ic,mz]=mzscan(obj, mztarget, varargin)
    % Get ion counts for given mztarget across elution times
    % Also return weighted mean m/z for each elution time
      defaults=struct('mztol',0.01);
      args=processargs(defaults,varargin);
      % Build a table of elution vs IC for this M/Z
      ic=zeros(size(obj.peaks,1),1);
      mz=nan(size(obj.peaks,1),1);
      for i=1:length(obj.peaks)
        sel=abs(obj.peaks{i}(:,1)-mztarget)<=args.mztol;
        if any(sel)
          ic(i)=sum(obj.peaks{i}(sel,2));
          mz(i)=sum((obj.peaks{i}(sel,1).*obj.peaks{i}(sel,2)))/ic(i);
        end
      end
    end
    
    function res=findcompound(obj, mztarget, varargin)
    % Find compound directly (not using clusters)
    % Returns struct containing possible peaks (integrated over peak in both time and m/z) at different elutions times
    % Only peaks with ion count >= peakratio* maximum peak ion count are returned
    % If sdf provided, only used to construct a name for the compound
      defaults=struct('name','','sdf','','elutetime',[],'mztol',0.01,'timetol',30,'debug',false,'ignoreelutetimes',[],'peakratio',0.05);
      args=processargs(defaults,varargin);
      if isempty(args.name)
        if ~isempty(args.sdf)
          args.name=sprintf('%s.%s',args.sdf.BATCH_PLATE,args.sdf.BATCH_WELL);
        else
          args.name=sprintf('M/Z=%f',mztarget);
        end
      end
      
      % Build a table of elution vs IC for this M/Z
      [ic,mz]=obj.mzscan(mztarget,'mztol',args.mztol);

      % Handle any time limitations given on command line
      if ~isempty(args.elutetime)
        sel=abs(obj.time-args.elutetime)<=args.timetol;
      elseif ~isempty(args.ignoreelutetimes)
        sel=~(min(abs(obj.time-args.ignoreelutetimes(:)'),[],2)<args.timetol);
      else
        sel=true(size(ic));
      end

      time=obj.time(sel);
      ic=ic(sel);
      mz=mz(sel);
      
      % Find peaks
      maxtime=[]; maxic=[]; maxmz=[];
      if ~isempty(ic) && length(ic)>1
        [p,pfwhh]=mspeaks(time,ic,'OVERSEGMENTATIONFILTER',args.timetol);
        if ~isempty(p)
          for i=1:size(p,1)
            sel=time>=pfwhh(i,1) & time<=pfwhh(i,2);
            maxic(i)=sum(ic(sel));
            maxmz(i)=nansum(mz(sel).*ic(sel))/maxic(i);
            maxtime(i)=p(i,1);
          end
          [maxic,ord]=sort(maxic,'desc');
          maxmz=maxmz(ord);
          maxtime=maxtime(ord);
          keep=maxic>=max(maxic)*args.peakratio;
          maxic=maxic(keep);
          maxmz=maxmz(keep);
          maxtime=maxtime(keep);
        end
      end
      
      desc=sprintf('m/z=%.4f+=%.3f',mztarget,args.mztol);
      if ~isempty(args.elutetime)
        desc=[desc,sprintf(' T=%.1f+=%.1f',args.elutetime,args.timetol)];
      elseif ~isempty(args.ignoreelutetimes)
        desc=[desc,sprintf(' T!=[%s]+=%.1f',sprintf('%.1f,',args.ignoreelutetimes),args.timetol)];
      end

      res=struct('mztarget',mztarget,'desc',desc,'findargs',args,'mz',maxmz,'time',maxtime,'ic',maxic);
      if isempty(maxic)
        if args.debug
          fprintf('No match to %s\n',desc);
        end
      elseif length(maxmz)==1
        if args.debug
          fprintf('Unique peak for %s:  T=%.0f IC=%.0f\n',desc, maxtime(1),maxic(1));
        end
      else
        if args.debug
          fprintf('Have multiple peaks for %s with peak ratio=%f\n', desc, maxic(2)/maxic(1));
        end
      end
    end

    function listidents(obj)
    %  List all the identified compounds
      for i=1:length(obj.ident)
        id=obj.ident(i);
        fprintf('%2d %-12.12s M/Z=%8.4f\n',i, id.name,id.mztarget);
        for k=1:length(id.peaks)
          pk=id.peaks(k);
          fprintf('   %12.12s     %8.4f  T=[%6.1f,%6.1f,%6.1f]  IC=%8.0f\n',  '',pk.mzobs, pk.elutions, pk.ictotal);
        end
      end
    end

    function plotident(obj,i)
    % Make plot for i-th identification entry
      id=obj.ident(i);
      if isempty(id.peaks)
        fprintf('No peaks for ident(%d)\n',i);
        return;
      end
      
      setfig(sprintf(' %.4f',id.mztarget));clf;

      subplot(211);
      pk=id.peaks(1);
      semilogy(obj.time,id.ioncount);
      hold on;
      ax=axis;
      tmargin=ceil((pk.elutions(3)-pk.elutions(1))/4);
      ax(1)=obj.time(max(1,pk.elutions(1)-tmargin));
      ax(2)=obj.time(min(length(obj.peaks),pk.elutions(3)+tmargin));
      axis(ax);
      plot(ax(1:2),pk.icthresh*[1,1],':r');
      plot(obj.time(pk.elutions(1))*[1,1],ax(3:4),':r');
      plot(obj.time(pk.elutions(2))*[1,1],ax(3:4),':g');
      plot(obj.time(pk.elutions(3))*[1,1],ax(3:4),':r');
      xlabel('Time');
      ylabel('Ion count');
      title(sprintf('M/Z %.4f (err=%.4f) IC=%.0f',id.mztarget,pk.mzobs-id.mztarget,pk.ictotal));

      subplot(212);
      stem(obj.peaks{pk.elutions(2)}(:,1),obj.peaks{pk.elutions(2)}(:,2));
      xlabel('m/z');
      ylabel('Ion count');
      title(sprintf('T(%d)=%f',pk.elutions(2),obj.time(pk.elutions(2))));
    end


    function plotmz(obj,time)
    % Plot m/z vs ioncount for trace nearest given time
      if isempty(obj.resamp)
        obj.resample();
      end
      [~,closest]=min(abs(time-obj.time));
      ti=sprintf('m/z @ T=%f',obj.time(closest));
      setfig(ti);clf;
      plot(obj.resamp.mz,obj.resamp.y(:,closest));
      title(ti)
      xlabel('m/z')
      ylabel('Relative Intensity')
      hold on;
      ax=axis;
      for i=1:length(obj.ident)
        id=obj.ident(i);
        plot(id.mztarget*[1,1],[0,ax(4)*1.1],':r')
      end
      % List major peaks
      pks=obj.peaks{closest};
      [~,ord]=sort(pks(:,2),'desc');
      for ii=1:10
        i=ord(ii);
        fprintf('M/Z=%.4f, IC=%6.0f\n', pks(i,:));
      end
      
    end
    
    function ic=TIC(obj)
      ic=cellfun(@(z) sum(z(:,2)),obj.peaks);
    end
    
    function plotTIC(obj)
    % Plot TIC showing best traces
      setfig('TIC');clf;
      plot(obj.time,obj.TIC());
      title('Total Ion Chromatogram (TIC)')
      xlabel('Retention Time')
      ylabel('Relative Intensity')
      hold on;
      ax=axis;
      for i=1:length(obj.ident)
        id=obj.ident(i);
        for k=1:length(id.peaks)
          pk=id.peaks(k);
          if k==1
            col='g';  % Green for main peak
          else
            col='r';   % Red for others
          end
          plot(obj.time(pk.elutions(2))*[1,1],[0,ax(4)*1.1],[':',col])
        end
      end
    end

    function map=timealign(obj,o2,step)
    % Find mapping between elution times
      t1=obj.TIC();
      t2=o2.TIC();
      if nargin<3
        step=10;
        window=500;
      end
      map=[];
      for i=1:step:length(t1)-window
        ind=finddelay(t1(i:i+window-1).*hamming(window),t2(i:i+window-1).*hamming(window),window/10);
        map(end+1,:)=[i+window/2,i+window/2+ind];
      end
      setfig(sprintf('timealign %s-%s',obj.name,o2.name));
      plot(map(:,1),map(:,2)-map(:,1),'.-');
      xlabel(obj.name);
      ylabel(o2.name);
    end
    
    
    function x=compare(obj,o2)
    % Compare identified peaks with another MassSpec object
      [c,ia,ib]=intersect({obj.ident.name},{o2.ident.name});

      fprintf('Matched %d of %d+%d identified compounds\n', length(c), length(obj.ident), length(o2.ident));
      fprintf('                         %20.20s                      %20.20s\n', obj.name, o2.name);
      %fprintf('Name           M/Z Exp.   dM/Z Elution Ion Cnt        dM/Z Elution Ion Cnt  IC 2/1\n');
      elutions=[];
      x=table();
      for i=1:length(c)
        i1=obj.ident(ia(i));
        i2=o2.ident(ib(i));
        if i1.mztarget~=i2.mztarget
          fprintf('Warning: %s has different M/Z targets (%.3f, %.3f)\n', i1.mztarget, i2.mztarget);
        end

        peakdelta=1e10;
        pk1=1;pk2=1;
        for k1=1:length(i1.peaks)
          t1=obj.time(i1.peaks(k1).elutions(2));
          for k2=1:length(i2.peaks)
            t2=o2.time(i2.peaks(k2).elutions(2));
            if abs(t1-t2)<abs(peakdelta)
              pk1=k1;
              pk2=k2;
              peakdelta=t2-t1;
            end
          end
        end
        if peakdelta==1e10
          continue;
        end
        p1=i1.peaks(pk1);
        p2=i2.peaks(pk2);
        row=table();
        row.name=i1.name;
        row.mztarget=i1.mztarget;
        row.id1=sprintf('%4d.%-3d',ia(i),pk1);
        row.mzdelta1=round(p1.mzobs-i1.mztarget,3);
        row.elution1=obj.time(p1.elutions(2));
        row.ic1=round(p1.ictotal,0);

        row.id2=sprintf('%4d.%-2d',ib(i),pk2);
        row.mzdelta2=round(p2.mzobs-i1.mztarget,3);
        row.elution2=o2.time(p2.elutions(2));
        row.ic2=round(p2.ictotal,0);
        row.tdelta=round(peakdelta,2);
        row.flags='  ';
        if abs(p1.mzobs-p2.mzobs)>0.005
          row.flags(1)='M';
        end
        if abs(peakdelta)>10
          row.flags(2)='T';
        end
        %fprintf('%-12.12s   %7.3f  %6.3f %7.2f %7.0f      %6.3f %7.2f %7.0f %6.1f%%\n', c{i}, i1.mztarget, p1.mzobs-i1.mztarget, p1.elutions(2), p1.ictotal, i2.mzobs-i1.mztarget, p2.elutions(2), p2.ictotal, p2.ictotal/p1.ictotal*scaling*100);
        x=vertcat(x,row);
      end
      scaling=median(x.ic1./x.ic2);   % Scale factor to match ion counts
      x.icratio=round(100*x.ic2./x.ic1*scaling,1);
      disp(x)
      ti=sprintf('%s vs %s',obj.name, o2.name);

      setfig(ti);clf;
      subplot(221);
      loglog(x.ic1, x.ic2,'o');
      xlabel(sprintf('Ion Count %s',obj.name));
      ylabel(sprintf('Ion Count %s',o2.name));
      ax=axis;
      hold on;
      plot(ax(1:2),ax(1:2)/scaling,':r');
      title(ti);
      subplot(222);
      plot(x.mztarget, x.mzdelta1,'o');
      hold on;
      plot(x.mztarget, x.mzdelta2,'x');
      xlabel('M/Z Exact')
      ylabel('M/Z Diff')
      subplot(223);
      plot(x.elution1,x.elution2,'o');
      xlabel('T1 (sec)');
      ylabel('T2 (sec)')
    end
    
  end % methods
end % classdef