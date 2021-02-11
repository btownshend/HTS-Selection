classdef Chem < handle
  properties
  end
  
  methods(Static)
    function fs=formula2struct(f)
    % Convert to a struct
      assert(ischar(f));
      fs=struct();
      upper=find(f>='A' & f<='Z');
      upper(end+1)=length(f)+1;
      for i=1:length(upper)-1
        pos=upper(i);
        if pos<length(f) & f(pos+1)>='a' & f(pos+1)<='z'
          sym=f(pos:pos+1);
          npos=pos+2;
        else
          sym=f(pos);
          npos=pos+1;
        end
        if npos<upper(i+1)
          fs.(sym)=str2num(f(npos:upper(i+1)-1));
        else
          fs.(sym)=1;
        end
      end
      f=fs;
    end

    function f=struct2formula(fs)
      assert(isstruct(fs));
      fn=fieldnames(fs);
      f='';
      for i=1:length(fn)
        if fs.(fn{i})==1
          f=sprintf('%s%s',f,fn{i});
        else
          f=sprintf('%s%s%d',f,fn{i},fs.(fn{i}));
        end
      end
    end
    
    function c=getisotopes(f,varargin)
    % Get list of isotopes for given formula (a struct), returning a struct of name,mass,abundance (relative to original)
      defaults=struct('minabundance',1e-4,'fine',false,'plot',false);
      args=processargs(defaults,varargin);

      if ischar(f)
        f=Chem.formula2struct(f);
      end
      
      isoabund={1,'H',99.9885,1.007825   % n, element, %abund, mass
                2,'H',0.0115,2.0140
                12,'C',98.93,12
                13,'C',1.07,13.00335
                14,'N',99.636,14.00307
                15,'N',0.364,15.00011
                16,'O',99.757,15.99491
                17,'O',0.038,16.99913
                18,'O',0.205,17.99916
                19,'F',100.0,18.99840
                28,'Si',92.223,27.97693
                29,'Si',4.685,28.97649
                30,'Si',3.092,29.97376
                31,'P',100,30.9737619986
                32,'S',94.99,31.97207
                33,'S',0.75,32.97146
                34,'S',4.25,33.96786
                35,'Cl',75.76,34.96885
                37,'Cl',24.24,36.96590
                79,'Br',50.69,78.9183
                81,'Br',49.31,80.9163
               };
      isoabund=cell2struct(isoabund,{'n','element','abund','mass'},2);
      for i=1:length(isoabund)
        isoabund(i).abund=isoabund(i).abund/100;
        if isoabund(i).abund>0.5
          isoabund(i).name=isoabund(i).element;
        else
          isoabund(i).name=sprintf('[%d]%s',isoabund(i).n,isoabund(i).element);
        end
      end

      atoms=fieldnames(f);
      if length(atoms)>1
        % Multiple atoms, handle recursively
        c1=Chem.getisotopes(struct(atoms{1},f.(atoms{1})),'minabundance',args.minabundance);
        c2=Chem.getisotopes(rmfield(f,atoms{1}),'minabundance',args.minabundance);
        c=[];
        for i=1:length(c1)
          for j=1:length(c2)
            abund=c1(i).abundance*c2(j).abundance;
            if abund > args.minabundance
              c=[c,struct('name',[c1(i).name,c2(j).name],'mass',c1(i).mass+c2(j).mass,'abundance',abund)];
            end
          end
        end
      else
        % Only 1 atom type
        a=isoabund(strcmp({isoabund.element},atoms{1}));   
        if isempty(a)
          error('Bad element: %s', atoms{1});
        end
        assert(abs(sum([a.abund])-1)<= .001);
        
        n=f.(atoms{1});
        v=zeros(length(a),1);
        v(1)=n;
        c=[];
        while true
          abund=factorial(n);
          for i=1:length(a)
            abund=abund*(a(i).abund^v(i))/factorial(v(i));
          end
          %fprintf('v=[%s], a=%f\n', sprintf('%d,',v),abund);
          if abund>args.minabundance
            name='';
            for k=1:length(a)
              if v(k)>1
                name=sprintf('%s%s%d',name,a(k).name,v(k));
              elseif v(k)>0
                name=sprintf('%s%s',name,a(k).name);
              end
            end
            c=[c,struct('name',name,'mass',dot([a.mass],v),'abundance',abund)];
          else
            v(end)=9999;
          end
          vpos=length(v);
          while true
            v(vpos)=v(vpos)+1;
            v(1)=n-sum(v(2:end));
            if v(1)>=0
              break;
            end
            v(vpos)=0;
            vpos=vpos-1;
            if vpos==1
              break;
            end
          end
          if vpos==1
            break;
          end
        end
      end
      assert(abs(sum([c.abundance])-1)<.02);

      % If this is a coarse pattern, merge peaks for same integer difference
      if ~args.fine
        % Merge isotopes that are too close to resolve
        [~,ord]=sort([c.mass]);
        c=c(ord);
        i=1;
        while i<length(c)
          if c(i+1).mass-c(i).mass < 0.01
            abund=[c([i,i+1]).abundance];
            total=sum(abund);
            meanmass=(c(i).mass*sum(c(i).abundance)+c(i+1).mass*sum(c(i+1).abundance))/total;
            if sum(c(i+1).abundance)>sum(c(i).abundance)
              c=c([1:i-1,i+1:end]);
            else
              c=c([1:i,i+2:end]);
            end
            c(i).abundance=total;
            c(i).mass=meanmass;
          else
            i=i+1;
          end
        end
      end
      
      [~,ord]=sort([c.abundance],'desc');
      c=c(ord);
      
      if args.plot
        setfig('getisotopes');clf;
        stem([c.mass],[c.abundance]);
        set(gca,'YScale','log');
        xlabel('Mass');
        ylabel('Abundance');
      end
      if nargout<1
        struct2table(c)
      end
    end
    
    function c=dbsavecompound(name,formula)
    % Save compound in database along with isotopes and return PK
    %mysql('use hts');
      c=mysql(sprintf('select compound from compounds where id=''%s''',name));
      if ~isempty(c)
        fprintf('Database already contains %s\n', name);
        return;
      end
      mysql('start transaction');
      if isstruct(formula)
        formula=Chem.struct2formula(formula);
      end
      f=mysql(sprintf('select formula_pk from formulas where formula=''%s''',formula));
      if isempty(f)
        mysql(sprintf('insert into formulas(formula) values(''%s'')',formula));
        f=mysql(sprintf('select formula_pk from formulas where formula=''%s''',formula));
        iso=Chem.getisotopes(formula);
        vals='';
        for i=1:length(iso)
          vals=sprintf('%s\n(%d,''%s'',%f,%f),',vals,f,iso(i).name,iso(i).mass,iso(i).abundance);
        end
        vals=vals(1:end-1);  % Remove trailing comma
        mysql(sprintf('insert into isotopes(formula_pk,formula,mass,abundance) VALUES\n%s',vals));
        fprintf('Added %d isotopes for %s\n', length(iso), formula);
      end
      mysql(sprintf('insert into compounds(id,formula_pk) values(''%s'',%d)',name,f));
      c=mysql(sprintf('select compound from compounds where id=''%s''',name));
      mysql('commit');
      fprintf('Created compound %d for %s (%s)\n', c, name, formula);
    end
    
  end
end % Chem