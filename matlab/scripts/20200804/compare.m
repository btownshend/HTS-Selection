% Check dilution series
msdir='../../../data/MassSpec/';
setup={
'Blank','DMSO',0,1,
'Blank','Water',0,1,
'Blank','ACN',0,1,
'CDIV','DMSO',200,1,
'CDIV','DMSO',100,1,
'CDIV','DMSO',50,1,
'CDIV','DMSO',200,2,
'CDIV','DMSO',100,2,
'CDIV','DMSO',50,2,
'CDIV','DMSO',200,3,
'CDIV','DMSO',100,3,
'CDIV','DMSO',50,3,
'CDIV','Water',200,1,
'CDIV','Water',100,1,
'CDIV','Water',50,1,
'CDIV','Water',200,2,
'CDIV','Water',100,2,
'CDIV','Water',50,2,
'CDIV','Water',200,3,
'CDIV','Water',100,3,
'CDIV','Water',50,3,
'CDIV','ACN',200,1,
'CDIV','ACN',100,1,
'CDIV','ACN',50,1,
'CDIV','ACN',200,2,
'CDIV','ACN',100,2,
'CDIV','ACN',50,2,
'CDIV','ACN',200,3,
'CDIV','ACN',100,3,
'CDIV','ACN',50,3,
'Blank','DMSO',0,2,
    }
folder=[msdir,'20200803'];
if ~exist('dilms','var')
  dilms={};
end
for i=1:size(setup,1);
  if i> length(dilms)
    s=setup(i,:);
    name=s{1}; solvent=s{2}; conc=s{3}; repl=s{4};
    if conc==0
      path=sprintf('%s/%s_%s_%d.mzXML',folder,name,solvent,repl);
    else
      path=sprintf('%s/%s_%s_%.0fnM_%d.mzXML',folder,name,solvent,conc,repl);
    end
    dilms{i}=MassSpec(path);
    dilms{i}.setLoad(10e-6*conc*1e-9);
    dilms{i}.name=sprintf('%s-%s@%.0f.%d',name,solvent,conc,repl);
    % Prune out some data
    %dilms{i}.filter([200,2900],[130,530]);  % NOTE: this is using the localtimes to filter
  end
end

mz=compounds.mass+compounds.ADDUCTS(1).mass;
mz=mz(:);
ref=13;
comp=MassSpec.comparepeaks(dilms,'mztol',.01,'minic',5000,'mz',mz,'ref',ref);
lcomp=log(comp(:,2:end));
lcomp(comp(:,2:end)==0)=nan;
cc=corrcoef(lcomp,'rows','pairwise');
data=cc; data(end+1,:)=nan; data(:,end+1)=nan;
setfig('correlation');clf;
pcolor(data);
set(gca,'XTick',(1:length(dilms))+0.5);
set(gca,'XTickLabel',cellfun(@(z) z.name, dilms,'Unif',false));
set(gca,'XTickLabelRotation',90);
set(gca,'YTick',(1:length(dilms))+0.5);
set(gca,'YTickLabel',cellfun(@(z) z.name, dilms,'Unif',false));
colorbar;

setfig('Scaling');clf;
cm=nanmedian(comp(:,find(strcmp(setup(:,1),'CDIV'))+1),2);
r=comp(:,2:end);
ratio=[];
for i=1:size(r,2)
  r(:,i)=r(:,i)./cm;
  ratio(i)=nanmedian(r(:,i));
end
conc=[];
for i=1:size(setup,1)
  conc(i)=setup{i,3};
end
conc=conc';
name=setup(:,2);
uname=unique(name);
h=[];
for i=1:length(uname)
  sel=strcmp(name,uname{i});
  h(i)=plot(conc(sel),ratio(sel),'o');
  hold on;
  uconc=unique(conc(sel))';
  meanr=[];
  for j=1:length(uconc)
    meanr(j)=mean(ratio(sel & conc==uconc(j)));
  end
  plot(uconc,meanr,'-','Color',get(h(i),'Color'));
  % Also plot straight line
  avg=sum(meanr)/sum(uconc);
  plot(uconc,avg*uconc,':','Color',get(h(i),'Color'));
  keyboard;
end
legend(h,uname);
xlabel('Concentration');
ylabel('Relative Ion Counts');


