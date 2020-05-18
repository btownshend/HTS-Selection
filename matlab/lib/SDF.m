% Class for dealing with SDF files
classdef SDF < handle
    properties
      sdf;
    end

    properties(Constant)
      valence=struct('C',4,'O',2,'H',1,'N',[3,5],'F',1,'Cl',1,'Br',1,'S',[2,4,6],'P',[3,5],'B',3);
    end
    
    methods
      function obj=SDF()
        ;
      end

      function n=getname(obj,i,withnum)
        n=sprintf('%s%s',strrep(obj.sdf(i).BATCH_PLATE,'CDIV0',''),obj.sdf(i).BATCH_WELL);
        if n(1)=='0'
          n=n(2:end);
        end
        if nargin<3 || withnum
          n=sprintf('%s(%d)',n,i);
        end
      end
      
      function f=getformula(obj,i)
        f=[];
        F=obj.sdf(i).Formula;
        fn=fieldnames(F);
        for i=1:length(fn)
          f=[f,fn{i}];
          if F.(fn{i})~=1
            f=[f,num2str(F.(fn{i}))];
          end
        end
      end
      
      function sel=find(obj,plate,row,col)
        if iscell(plate)
          sel=obj.find(plate{1});
          for i=1:length(plate)
            sel=sel|obj.find(plate{i});
          end
          return;
        end
        if nargin==2 && ischar(plate)
          ind=find(plate>='A' & plate<='H');
          if length(ind)~=1
            error('Unable to parse %s (expected NNNNANN)\n', plate);
          end
          sel=obj.find(str2num(plate(1:ind-1)),plate(ind),str2num(plate(ind+1:end)));
          return;
        end
        
        sel=ones(size(obj.sdf));
        if ~isempty(plate)
          sel=sel&arrayfun(@(z) ismember(str2num(z.BATCH_PLATE(5:end)),plate),obj.sdf);
        end
        if nargin>=3 && ~isempty(row)
          sel=sel&arrayfun(@(z) ismember(z.BATCH_WELL(1),row),obj.sdf);
        end
        if nargin>=4 && ~isempty(col)
          sel=sel&arrayfun(@(z) ismember(str2num(z.BATCH_WELL(2:end)),col),obj.sdf);
        end
      end
      
      function newsdf=filter(obj,sel)
        newsdf=SDF();
        newsdf.sdf=obj.sdf(sel);
      end
      
      function read(obj,file,maxrecords)
      % Load an SDF file into matlab
        obj.sdf=[];
        fd=fopen(file,'r');
        if fd<0
          error('Unable to open %s for reading\n', file);
        end
        fprintf('Loading %s...',file);
        nread=0;
        while true
          if mod(nread,100)==0
            fprintf('%d...',length(obj.sdf));
          end
          header=textscan(fd,'%s',3,'Delimiter','\n');
          header=header{1};
          [~,enum]=ferror(fd);
          if enum~=0
            break;
          end
          line=fgetl(fd);
          s=textscan(line,'%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%6s',1);
          natoms=s{1};
          nbonds=s{2};
          
          atoms=[];
          for i=1:natoms
            line=fgetl(fd);
            s=textscan(line,'%10f%10f%10f %3s%3d%3d%3d',1);
            atoms=[atoms,struct('x',s{1},'y',s{2},'z',s{3},'atom',s{4},'massdiff',s{5},'charge',s{6})];
          end
          bonds=[];
          for i=1:nbonds
            line=fgetl(fd);
            s=textscan(line,'%3d%3d%3d%3d%3d%3d%3d',1);
            bonds=[bonds,struct('atoms',[s{1},s{2}],'type',s{3},'stereo',s{4})];
          end
          line=fgetl(fd);
          while ~strcmp(line,'M  END')
            % Skip charge lines
            %fprintf('Skip ''%s''\n',line);
            line=fgetl(fd);
          end
          assert(strcmp(line,'M  END'));
          line=fgetl(fd);
          while (isempty(line))
            line=fgetl(fd);
          end
          obj.sdf(end+1).header=header;
          obj.sdf(end).atoms=atoms;
          obj.sdf(end).bonds=bonds;
          while line(1)=='>'
            dataname=textscan(line(3:end),'<%[^>]>%s',1);
            dataval='';
            line=fgetl(fd);
            while line(1)~='>' && ~strcmp(line,'$$$$')
              dataval=[dataval,line];
              line=fgetl(fd);
              while (isempty(line))
                line=fgetl(fd);
              end
            end
            nval=str2num(dataval);
            if ~isempty(nval)
              dataval=nval;
            end
            obj.sdf(end).(dataname{1}{1})=dataval;
          end
          assert(strcmp(line,'$$$$'));
          nread=nread+1;
          if nargin>1 && nread>=maxrecords
            break;
          end
        end
        fclose(fd);
        fprintf('%d done\n',nread);
      end

      
      function getformulae(obj,force)
      % Get formulae for each entry
      % Use same format as 'F' in isotopicdist
        fprintf('Computing formulae and masses...');
        if nargin<2
          force=false;
        end
        for i=1:length(obj.sdf)
          if ~force && isfield(obj.sdf(i),'Formula') && ~isempty(obj.sdf(i).Formula)
            continue;
          end
          if mod(i,100)==0
            fprintf('%d...',i);
          end
          % Compute the explicit structure
          f=struct('C',0,'H',0);
          atoms=obj.sdf(i).atoms; bonds=obj.sdf(i).bonds;
          [anames,ia,ic]=unique({atoms.atom});
          for j=1:length(anames)
            f.(anames{j})=sum(ic==j);
          end
          % Add missing Hydrogens
          % "The implicit hydrogen count is determined by summing the bond orders of the bonds connected to the atom. If that sum is equal to a known valence for the element or is greater than any known valence then the implicit hydrogen count is 0. Otherwise the implicit hydrogen count is the difference between that sum and the next highest known valence."
          missingValence=0;
          for j=1:length(atoms)
            val=sort(obj.valence.(atoms(j).atom));
            sel=arrayfun(@(z) any(z.atoms==j),bonds);
            nbonds=sum([bonds(sel).type]);
            if ~ismember(nbonds,val)
              nextVal=val(find(val>nbonds,1));
              if ~isempty(nextVal)
                missingValence=missingValence+nextVal-nbonds;
              end
            end
          end
          f.H=f.H+missingValence;
          while true
            [~,info]=isotopicdist(f);
            fn=fieldnames(info);
            for j=1:length(fn)
              obj.sdf(i).(fn{j})=info.(fn{j});
            end
            expectedAvgMass=obj.sdf(i).BATCH_MW;
            if ~isempty(obj.sdf(i).Addend_Display) && ~strcmp(obj.sdf(i).Addend_Display,'Hydrate')
              expectedAvgMass=expectedAvgMass-obj.sdf(i).Addend_MolWeight*obj.sdf(i).BatAdd_Equivalents;
            end
            extraHs = expectedAvgMass-obj.sdf(i).ObservedAverageMass;
            if abs(extraHs)>2.5 || extraHs+f.H < -0.5
              fprintf('Mass difference for structure %s (%s) - expected mw=%.1f, computed=%.1f - off by %.1f\n',obj.getname(i),obj.getformula(i),expectedAvgMass, obj.sdf(i).ObservedAverageMass, extraHs);
              keyboard;
            elseif extraHs>0.5
                fprintf('Mass difference for structure %s (%s) - expected=%.1f, computed=%.1f - adding %d H\n',obj.getname(i),obj.getformula(i),expectedAvgMass,obj.sdf(i).ObservedAverageMass,round(extraHs));
                f.H=f.H+round(extraHs);
                continue;
            elseif extraHs<-0.5
                fprintf('Mass difference for structure %s (%s) - expected=%.1f, computed=%.1f - dropping %d H\n',obj.getname(i),obj.getformula(i),expectedAvgMass,obj.sdf(i).ObservedAverageMass,round(-extraHs));
                f.H=f.H+round(extraHs);
                continue;
            end
            break;
          end
        end
        fprintf('done\n');
      end
      
      function csvwrite(obj,file)
        fd=fopen(file,'w');
        fprintf(fd,'Plate,Row,Col,ID,Formula,MW,SaltMW\n');
        for i=1:length(obj.sdf)
          s=obj.sdf(i);
          fprintf(fd,'%d,%s,%d,%s,%s,%.4f,%.4f\n', str2num(s.BATCH_PLATE(5:end)),s.BATCH_WELL(1),str2num(s.BATCH_WELL(2:end)),s.compound_Corp_Reg_Number,obj.getformula(i),s.MonoisotopicMass,s.saltMass);
        end
        fclose(fd);
      end

      function plot(obj,sel,pos,align,aspect)
        if islogical(sel)
          sel=find(sel);
        end
        if length(sel)>1
          % Multiple structures -- layout in grid
          %clf;
          s=obj.sdf(sel);
          maxsize=[0,0,0];
          for i=1:length(s)
            [low,high]=bounds([s(i).atoms.x;s(i).atoms.y;s(i).atoms.z]');
            maxsize=max(maxsize,high-low);
          end
          nz=1; % ceil(length(sel)^(1/3));
          if nargin<5
            aspect=1;
          end
          nx=ceil(sqrt(aspect*length(sel)/nz));
          ny=ceil(length(sel)/nx/nz);
          i=1;
          shift=maxsize+[0.5,1.5,2];
          for iy=ny:-1:1
            for ix=1:nx
              for iz=1:nz
                %subplot(nr,nc,i);
                if i<=length(sel)
                  obj.plot(sel(i),shift.*[ix,iy,iz],'lowermid');
                end
                i=i+1;
              end
            end
          end
          
          return;
        end
        
        if isempty(sel)
          fprintf('Nothing to plot\n');
          return;
        end
        
        if nargin<4 || isempty(align)
          align='center';
        end
        if nargin<3 || isempty(pos)
          pos=[0,0,0];
        end
        s=obj.sdf(sel);
        [low,high]=bounds([s.atoms.x;s.atoms.y;s.atoms.z]');
        if strcmp(align,'lowerleft')
          pos=pos-low;
        elseif strcmp(align,'lowermid')
          pos(2:3)=pos(2:3)-low(2:3);
          pos(1)=pos(1)-(high(1)+low(1))/2;
        elseif strcmp(align,'center');
          pos=pos-(high+low)/2;
        else
          error('Unsupported alignment: %s\n', align);
        end

        b=s.bonds;
        for i=1:length(b)
          a=s.atoms(b(i).atoms);
          v=[[a.x];[a.y];[a.z]]';
          delta=v(2,:)-v(1,:);
          if a(1).atom=='C'
            vp=v(1,:)+pos;
          else
            vp=v(1,:)+pos+delta*0.2;
          end
          if a(2).atom=='C'
            vp(2,:)=v(1,:)+pos+delta;
          else
            vp(2,:)=v(1,:)+pos+delta*0.8;
          end
          if b(i).type==1
            plot3(vp(:,1),vp(:,2),vp(:,3),'b-');
          else
            perp=delta([2,1,3]).*[1,-1,0];perp=perp/norm(perp)*0.07;
            assert(abs(dot(perp,delta))<1e-8);
            plot3(vp(:,1)-perp(1),vp(:,2)-perp(2),vp(:,3)-perp(3),'r-');
            plot3(vp(:,1)+perp(1),vp(:,2)+perp(2),vp(:,3)+perp(3),'r-');
          end
          hold on;
        end
        axis equal
        axis off
        view(0,90);
        %plot3([s.atoms.x],[s.atoms.y],[s.atoms.z],'.');
        noncarbon=s.atoms(~strcmp({s.atoms.atom},'C'));
        text([noncarbon.x]+pos(1),[noncarbon.y]+pos(2),[noncarbon.z]+pos(3),{noncarbon.atom},'HorizontalAlignment','center','VerticalAlignment','middle');
        text(mean([low(1),high(1)])+pos(1),pos(2)+low(2)-0.15,pos(3)+mean([low(3),high(3)]),obj.getname(sel,false),'HorizontalAlignment','center','VerticalAlignment','top','Color','m');
        %title(obj.getname(i));
      end
      
    end

  end