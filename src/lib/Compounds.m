% Data structure to hold information about compounds located in mass spec runs
classdef Compounds < handle
  properties
    names;   % names{i} - Name of compound i
    sdf;     % sdf{i} - SDF data for compound i
    mztarget; % mztarget(i) - target m/z for compound i 
    files;   % files{j} - Mass spec filename j
    moles;    % moles(j) - Moles of each compound loaded in run j
    group;   % group{j} - name of group this file belongs to
    contains;  % contains(i,j) is true if we expected to find compound i in file j
             % For compounds uniquely located in a particular run:
    mz;    % mz(i,j) contains the observed m/z peak on compound i, from file j
    time;    % time(i,j) contains the elution time on compound i, from file j
    ic;    % ic(i,j) contains the total ion count for compound i, from file j
    nisomers;   % nisomers(i,j) is the number of other compounds that are within MZFUZZ in this file
    nunique;   % nunique(i,j) is the number of other compounds that are within MZFUZZ in this file and with unknown elution time
    numhits;   % numhits(i,j) is the number of possible hits (need to have exactly 1 for the hit to be used)
  end
  
  properties(Constant)
    MZFUZZ=0.003;
    TIMEFUZZ=50;   % in seconds
  end
  
  methods
    function obj=Compounds()
    end

    function nindex=lookupName(obj,name, mztarget, sdf)
    % Find the index of compound by name, create if missing
      if ismember(name,obj.names)
        nindex=find(strcmp(name,obj.names));
      else
        obj.names{end+1,1}=name;
        nindex=length(obj.names);
        if nargin>=4
          obj.sdf{nindex,1}=sdf;
        else
          obj.sdf{nindex,1}=[];
        end
        if nargin>=3
          obj.mztarget(nindex,1)=mztarget;
        else
          obj.mztarget(nindex,1)=nan;
        end
        obj.mz(nindex,:)=nan;
        obj.time(nindex,:)=nan;
        obj.ic(nindex,:)=nan;
        obj.contains(nindex,:)=false;
      end
    end

    function findex=lookupMS(obj,ms)
    % Find the index of a file by name, create if missing
      if ismember(ms.name,obj.files)
        findex=find(strcmp(ms.name,obj.files));
      else
        obj.files{1,end+1}=ms.name;
        findex=length(obj.files);
        obj.mz(:,findex)=nan;
        obj.time(:,findex)=nan;
        obj.ic(:,findex)=nan;
        obj.contains(:,findex)=false;
        obj.moles(findex)=ms.moles;
      end
    end

    function addFromSDF(obj,ms,sdf,varargin)
    % Add the unique peaks for compounds in given SDF from a particular M/S run
    % Use prior analyses to figure out the expected elution time for each compound
    % or scan all elution times if the no prior data (keep only if a unique peak is determined)
      defaults=struct('debug',false,'group','');  
      args=processargs(defaults,varargin);

      fprintf('Adding data from %s\n',ms.name);
      findex=obj.lookupMS(ms);
      if ~isempty(args.group)
        obj.group{findex}=args.group;
      end
      % Add all these compounds and mark which ones this file contains
      nindex=[];
      for i=1:length(sdf.sdf)
        s=sdf.sdf(i);
        name=[s.BATCH_PLATE,'-',s.BATCH_WELL];
        mztarget=s.MonoisotopicMass+1;
        nindices(i)=obj.lookupName(name,mztarget,s);
        nindex=nindices(i);
        obj.contains(nindex,findex)=true;   % Mark it as expected to contain
      end
      % Attempt to locate each one uniquely
      for i=1:length(nindices)
        nindex=nindices(i);
        meantime=nanmean(obj.time(nindex,:));
        mztarget=obj.mztarget(nindex);
        % Check if the M/Z is unique over the compounds in this file
        samemz=nindices(find(abs(obj.mztarget(nindices)-mztarget)<=obj.MZFUZZ));   % Compounds with same M/Z
        nisomers=1;  nunique=1; ignoreelutetimes=[];
        for i=1:length(samemz)
          if samemz(i)~=nindex
            nisomers=nisomers+1;
            if ~any(isfinite(obj.time(samemz(i),:)))
              % This isomer has unknown elution time
              nunique=nunique+1;
            else
              ignoreelutetimes(end+1)=nanmean(obj.time(samemz(i),:));
            end
          end
        end
        obj.nisomers(nindex,findex)=nisomers;
        obj.nunique(nindex,findex)=nunique;
        
        if isfinite(meantime)
          id=ms.findcompound(mztarget,'sdf',s,'elutetime',meantime,'timetol',obj.TIMEFUZZ,'mztol',obj.MZFUZZ,'debug',args.debug);
        elseif nisomers==1
          id=ms.findcompound(mztarget,'sdf',s,'mztol',obj.MZFUZZ,'timetol',obj.TIMEFUZZ,'debug',args.debug);
        elseif nunique==1
          % Have multiple isomers but the other ones all have elution times already
          % Could look for addition peaks
          id=ms.findcompound(mztarget,'sdf',s,'mztol',obj.MZFUZZ,'timetol',obj.TIMEFUZZ,'debug',args.debug,'ignoreelutetimes',ignoreelutetimes);
        else
          if args.debug
            fprintf('Have %d indistinguishable isomers for %s in %s\n', nisomers, obj.names{nindex}, obj.files{findex});
          end
          continue;
        end
        obj.numhits(nindex,findex)=length(id.mz);
        if length(id.mz)==1
          obj.mz(nindex,findex)=id.mz(1);
          obj.time(nindex,findex)=id.time(1);
          obj.ic(nindex,findex)=id.ic(1);
        elseif length(id.mz)==0
          obj.ic(nindex,findex)=0;
        else
          ;
        end
      end
      fprintf('Located %d/%d compounds uniquely with %d not present in MS peaks\n', sum(isfinite(obj.mz(:,findex))), sum(obj.contains(:,findex)), sum(obj.ic(:,findex)==0));
    end
    
    function summary(obj)
    % Summarize data available
      fprintf('Contains %d files, %d compounds (%d with elution time)\n', length(obj.files), length(obj.names), sum(any(isfinite(obj.time'))));
      for i=1:length(obj.files)
        fprintf('%2d %-20.20s %3d/%3d compounds identified, %d missing\n', i, obj.files{i}, sum(isfinite(obj.mz(:,i))),sum(obj.contains(:,i)),sum(obj.nunique(:,i)==1 & obj.contains(:,i) & ~isfinite(obj.mz(:,i))));
      end
    end
    
    function x=report(obj)
    % Build table on data by compound
    % Each row is a single compound
    % Data by group
      ugroups=unique(obj.group,'sorted');
      x=[];
      [~,ord]=sort(obj.names);
      for ii=1:length(obj.names)
        i=ord(ii);
        x(ii).name=obj.names{i};
        x(ii).mztarget=obj.mztarget(i);
        
        for j=1:length(ugroups)
          files=find(strcmp(obj.group,ugroups{j})& obj.contains(i,:));
          x(ii).mzoffset(j)=nanmean(obj.mz(i,files))-obj.mztarget(i);
          x(ii).elution(j)=nanmean(obj.time(i,files));
          x(ii).ioncount(j)=nanmean(obj.ic(i,files));
          x(ii).nisomers(j)=nanmin(obj.nisomers(i,files));
          x(ii).nunique(j)=nanmin(obj.nunique(i,files));
          x(ii).numhits(j)=nanmin(obj.numhits(i,files));
          f='';
          for k=1:length(files)
            f=[f,obj.files{files(k)},','];
          end
          f=f(1:end-1);  % Remove trailing comma
          x(ii).files{j}=f;
        end
      end
      x=struct2table(x);
    end
    
    function pcolorplot(obj)
    % Heat map of ioncount in matrix of compounds * files
      data=obj.ic;
      % Normalize for xxth percentage
      for i=1:size(data,2)
        data(:,i)=data(:,i)/nanmean(data(:,i));
      end
      data(data<.01)=.01;
      data(data>10)=10;
      setfig('compounds');clf;
      data(end+1,:)=nan;
      data(:,end+1)=nan;
      pcolor(log10(data)');
      shading flat;
      colorbar;
      xlabel('Compound');
      ylabel('File');
      set(gca,'YTick',(1:length(obj.files))+0.5);
      set(gca,'YTickLabels',strrep(obj.files,'.mzXML',''));
      set(gca,'ticklabelinterpreter','none');
    end
    
    function scaling=getscaling(obj,f1,f2)
    % Get ioncounts in f2 relative to f1
      ratio=obj.ic(:,f2)./obj.ic(:,f1);
      scaling=nanmedian(ratio(ratio>0));
    end
    
    function checkmzoffset(obj)
    % Check whether mzoffset used when reading mass spec files should be changed
      fprintf('File          Additional Offset\n');
      ugroups=unique(obj.group);
      for j=1:length(ugroups)
        all=[];
        for i=1:length(obj.files)
          if strcmp(obj.group{i},ugroups{j})
            err=nanmedian(obj.mz(:,i)-obj.mztarget);
            fprintf('%-20.20s  %8.4f\n', obj.files{i}, err);
            all(end+1)=err;
          end
        end
        sel=strcmp(obj.group,ugroups{j});
        fprintf('%-20.20s  %8.4f over group\n', ugroups{j}, nanmedian(all));
      end
    end
    
    function plotcompare(obj,f1,f2)
    % Plot comparison of each compound the occurs in both f1 and f2
      ti=sprintf('%s vs %s',obj.files{f1}, obj.files{f2});
      setfig(ti);clf;
      subplot(221)
      plot(obj.mz(:,f1)-obj.mztarget,obj.mz(:,f2)-obj.mztarget,'o');
      hold on;
      ax=axis;
      plot(ax(1:2),obj.MZFUZZ*[1,1],'r:');
      plot(ax(1:2),-obj.MZFUZZ*[1,1],'r:');
      ax=axis;
      plot(obj.MZFUZZ*[1,1],ax(3:4),'r:');
      plot(-obj.MZFUZZ*[1,1],ax(3:4),'r:');
      xlabel(obj.files{f1},'Interpreter','none');
      ylabel(obj.files{f2},'Interpreter','none');
      title('m/z offset');

      subplot(222);
      plot(mean(obj.time(:,[f1,f2]),2),diff(obj.time(:,[f1,f2]),[],2),'o');
      hold on;
      ax=axis;
      plot(ax(1:2),obj.TIMEFUZZ*[1,1],'r:');
      plot(ax(1:2),-obj.TIMEFUZZ*[1,1],'r:');
      xlabel('Mean (s)');
      ylabel('Diff (s)');
      title('Elution Times');
      
      ratio=obj.ic(:,f2)./obj.ic(:,f1);
      scaling=nanmedian(ratio(ratio>0));
      ratio=ratio/scaling;
      ratio(ratio==0)=.001;
      ratio(obj.ic(:,f1)==0 & obj.ic(:,f2)~=0)=1000;

      subplot(223);
      x=obj.ic(:,f1); y=obj.ic(:,f2);
      x(x==0)=1;
      y(y==0)=1;
      loglog(x,y,'o');
      hold on;
      ax=axis;
      plot(ax(1:2),ax(1:2)*scaling,':');
      xlabel(obj.files{f1},'Interpreter','none');
      ylabel(obj.files{f2},'Interpreter','none');
      title('Ion Counts');
      
      subplot(224);
      %      semilogy(obj.mztarget,max(ratio,.01),'o');
      histogram(log10(ratio),20);
      xlabel('log10(IC2/IC1)');
      ylabel('N');
      title(sprintf('Ion Count Ratio (scaling=%.2f)',scaling));
      
      h=suptitle(sprintf('%s N=%d',ti,sum(all(isfinite(obj.ic(:,[f1,f2])),2))));
      set(h,'Interpreter','none');
    end
  end
end
