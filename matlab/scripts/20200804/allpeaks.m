% Check all DMSO peaks to figure out which ones are real
expts=1:length(dilms);
timerange=[1100,1120];

lib=setup(expts,1);
solvent=setup(expts,2);
conc=cell2mat(setup(expts,3));
repl=cell2mat(setup(expts,4));

if ~exist('pks','var')
  pks=MassSpec.comparepeaks(dilms(expts),'mztol',0.01,'minic',5000,'ref',4,'timerange',timerange);
end
ptotal=sum(pks(:,2:end),2);
[ptotal,ord]=sort(ptotal,'descend');
pks=pks(ord,:);


% Locate peaks that scale with concentration
ei=[];ed=[];
cr=cell2mat(setup(expts,3))';
cr=cr/mean(cr);
for i=1:size(pks,1)
  s=mean(pks(i,2:end));
  ei(i)=sum((pks(i,2:end)-s).^2);
  ed(i)=sum((pks(i,2:end)-s*cr).^2);
end
dep=ed<ei;   % Selector for concentration-dependent peaks
fprintf('Fraction concentration-dependent: %.2f\n', mean(dep));
dpks=pks(dep,:);
logdpks=log(dpks);
logdpks(dpks==0)=nan;
setfig('dpks'); clf;
loglog(dpks(1:min(end,1000),4+1),dpks(1:min(end,1000),5+1),'.');
xlabel(dilms{4}.name,'interpreter','none');
ylabel(dilms{5}.name,'interpreter','none');

% Plot fraction dependent on concentration
fdep=sum(dpks(:,2:end))./sum(pks(:,2:end));
setfig('FracDep');clf;
bar(fdep);
set(gca,'XTick',1:length(expts));
set(gca,'XTickLabel',cellfun(@(z) z.name, dilms(expts),'Unif',false));
set(gca,'XTickLabelRotation',45);
title('Fraction of ion count dependent on concentration');
ylabel('Fraction');

% Merge by concentration
uconc=unique(conc);

byconc=[];
for i=1:length(uconc)
  byconc(:,i)=nanmean(dpks(:,find(conc==uconc(i))+1),2);
end
setfig('byconc');clf;
plot(uconc,nanmedian(byconc));
xlabel('Conc (nM)');
ylabel('Mean Ion Count');

% Plot replicate 1 vs other replicates
setfig('Replicates');clf;
tiledlayout('flow');
for i=1:length(lib)
  if repl(i)==1
    others=find(strcmp(lib,lib{i}) & strcmp(solvent,solvent{i}) & conc==conc(i) & repl>repl(i));
    if ~isempty(others)
      nexttile;
      loglog(dpks(:,i+1),dpks(:,others+1),'.');
      p={};
      for kk=1:length(others)
        k=others(kk);
        cc=corrcoef(logdpks(:,k+1),logdpks(:,i+1),'Rows','pairwise');
        p{end+1}=sprintf('p=%.3f',cc(1,2));
      end
      xlabel('Repl 1');
      ylabel('Other replicates');
      title(dilms{i}.name);
      legend(p);
      axis equal;
    end
  end
end

% Plot solvents against each other
setfig('Solvents');clf;
tiledlayout('flow');
usolvent=unique(solvent);
sdpks=[];
for i=1:length(usolvent)
  sdpks(:,i)=sum(dpks(:,find(strcmp(solvent,usolvent{i}))+1),2);
end;
for i=1:length(usolvent)
  for j=i+1:length(usolvent)
    nexttile;
    loglog(max(sdpks(:,i),1e3),max(sdpks(:,j),1e3),'.');
    xlabel(usolvent{i});
    ylabel(usolvent{j});
    title(sprintf('Ratio=%.2f\n',sum(sdpks(:,j))/sum(sdpks(:,i))));
  end
end

setfig('Dropouts');
lgpks=sum(sdpks,2)>=1e5;
dropouts=mean(sdpks(lgpks,:)==0);
bar(dropouts);
set(gca,'XTick',1:length(usolvent));
set(gca,'XTickLabel',usolvent);
set(gca,'XTickLabelRotation',45);
ylabel('Fraction of large conc-dependent peaks missing');
xlabel('Solvent');
title('Dropouts');

% Figure out mappings
cnames=compounds.names;
mass=compounds.mass;

mz=[];
% C13 isotopes
c13=13.003355;
c12=12;
for i=1:length(mass)
  [m,n]=adducts(mass(i));
  mz=[mz;m',m'+c13-c12,m'+2*(c13-c12)];
  if i==1
    anames=[n;cellfun(@(z) [z,',C13'],n,'Unif',false);cellfun(@(z) [z,',C13*2'],n,'Unif',false)];
  end
end

mztol=0.007;
fd=fopen('matches.csv','w');
fprintf(fd,'M/Z\tIC\tOffset\tHits\n');
for i=1:min(5000,size(dpks,1))
  [dif,loc]=min(abs(dpks(i,1)-mz(:)));
  [k,j]=ind2sub(size(mz),loc);
  total=mz(loc);
  [mk,mj]=ind2sub(size(mz), find(abs(dpks(i,1)-mz(:))<=mztol));
  fprintf(fd,'%9.4f\t%6.2e\t%9.4f\t', dpks(i,1), sum(dpks(i,2:end)), dpks(i,1)-total);
  for mm=1:length(mk)
    if mm>1
      fprintf(fd,',');
    end
    fprintf(fd,'%s(%s)', cnames{mk(mm)},anames{mj(mm)});
  end
  fprintf(fd,'\n');
end
fclose(fd);

