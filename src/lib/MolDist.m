classdef MolDist < SDF
  properties
    dist;
  end
  
  methods
    function obj=MolDist()
      ;
    end
    
    function loadcsv(obj,file)
      fd=fopen(file,'r');
      line=fgetl(fd);
      labels=strsplit(line,',');
      labels=labels(2:end);
      fclose(fd);
      dist=csvread(file,1,1);  % Skip labels
      assert(length(labels)==size(dist,1));
      assert(length(labels)==size(dist,2));
      ord=[];
      names=arrayfun(@(z) obj.getname(z,false),1:length(obj.sdf),'UniformOutput',false);
      for i=1:length(labels)
        l=labels{i};
        if l(1)=='0'
          l=l(2:end);
        end
        ind=find(strcmp(l,names));
        if isempty(ind)
          error('Unable to find %s in SDF\n', l);
        end
        ord(i)=ind;
      end
      obj.dist=dist(ord,ord);   % Reorder
      assert(all(diag(obj.dist)==1));
    end

    function resort(obj)
    % Sort by plate, then col, then row
      names=arrayfun(@(z) obj.getname(z,false),1:length(obj.sdf),'UniformOutput',false);
      % Names are {p,pp,pp}rcc
      % Resorder to pppccr
      reord=cellfun(@(z) [sprintf('%03d',str2num(z(1:end-3))),z(end-1:end),z(end-3)], names, 'UniformOutput',false);
      [~,ord]=sort(reord);
      obj.dist=obj.dist(ord,ord);
      obj.sdf=obj.sdf(ord);
    end

    function setSDF(obj,s)
      obj.sdf=s.sdf;
    end
    
    function plotdist(obj)
      setfig('plotdist');clf;
      x=obj.dist;x(end+1,:)=nan;x(:,end+1)=nan;
      pcolor(x);
      shading flat;
      colorbar;
      ticks=1:80:length(obj.sdf);
      set(gca,'XTick',ticks);
      set(gca,'XTickLabel',arrayfun(@(z) obj.getname(z,false),ticks,'UniformOutput',false));
      set(gca,'XTickLabelRotation',45);
      set(gca,'YTick',ticks);
      set(gca,'YTickLabel',arrayfun(@(z) obj.getname(z,false),ticks,'UniformOutput',false));
      %set(gca,'YTickLabelRotation',45);
      colormap(hot);
      caxis([0,1]);
      title('Distance');
    end
    
    
    function x=closest(obj,names,varargin)
      defaults=struct('maxlist',20,'mindist',0.6);
      args=processargs(defaults,varargin);

      if ~iscell(names)
        names={names};
      end
      ind=[];
      for i=1:length(names)
        ind(i)=find(obj.find(names{i}));
      end
      d=obj.dist(ind,:);
      mdist=max(d,[],1);
      mdist(ind)=mdist(ind)+(length(ind):-1:1)/10000;  % Force initial ordering
      [mdist,ord]=sort(mdist,'desc');
      d=d(:,ord);
      onames=arrayfun(@(z) obj.getname(z,false),ord,'UniformOutput',false);
      fprintf('%6s ','');
      for i=1:length(names)
        fprintf('%6s ',names{i});
      end
      fprintf('\n');
      last=0;
      for i=1:min(args.maxlist,length(d))
        if mdist(i)<args.mindist
          break;
        end
        fprintf('%6s %s %s\n', onames{i}, sprintf('%6.2f ', d(:,i)),obj.getformula(ord(i)));
        last=i;
      end
      x=struct('colnames',{names},'dist',d(:,1:last)','index',ord(1:last)','rownames',{onames(1:last)});
      x.rownames=x.rownames';
    end
    
    function plotclosest(obj,x,ti)
    % Use results from closest and plot
      if isempty(ti)
        ti='plotclosest';
      end
      setfig(ti);clf;
      c=x.dist;
      c(end+1,:)=nan; c(:,end+1)=nan;
      pcolor(c);
      shading flat;
      caxis([0,1]);
      colorbar;
      set(gca,'XTick',(1:length(x.colnames))+0.5);
      set(gca,'XTickLabel',x.colnames);
      set(gca,'XTickLabelRotation',45);
      set(gca,'YTick',(1:length(x.rownames))+0.5);
      set(gca,'YTickLabel',x.rownames);
      colormap(hot);
      title(ti);
    end
  end 
end
