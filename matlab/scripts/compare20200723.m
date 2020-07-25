% Compare Massspec files from 7/23/20 test of various injections, elution times

datadir='../../data/';
matdir=[datadir,'matfiles/'];
resultsdir='../../results/';

msdir=[datadir,'MassSpec/'];


allfiles=dir([msdir,'20200723/*.mzXML']);


% Load data
if ~exist('testdata','var')
  testdata={};
  testmap={};
end

if length(testdata)>length(allfiles)
  fprintf('Removing last %d entries from testdata\n', length(testdata)-length(allfiles));
  testdata=testdata(1:length(allfiles));
end

for i=1:length(allfiles)
  path=[allfiles(i).folder,'/',allfiles(i).name];
  if i>length(testdata) || ~strcmp(testdata{i}.path,path)
      testdata{i}=MassSpec(path);
      if ~isempty(strfind(allfiles(i).name,'_10'))
        testdata{i}.setLoad(1000*1e-15);
      else
        testdata{i}.setLoad(500*1e-15);
      end
      % Prune out some data
      testdata{i}.filter([200,max(testdata{i}.time)-300],[130,530]); 
      testmap{i}=compounds.computeMap(testdata{i});
      pause(0.1);
  end
end

ref=mzdata{31};
refmap=maps(31);
for i=1:length(testdata)
  compounds.plotComposition(testdata{i},'map',testmap{i},'ref',ref,'refmap',refmap);
end

setfig('Injection Dependency');clf;
for i=1:2:length(testdata)
  subplot(2,3,(i+1)/2);
  testdata{i}.plotTIC(false);
  hold on;
  testdata{i+1}.plotTIC(false);
  l=legend(testdata{i}.name,testdata{i+1}.name,'Location','best');
  set(l,'Interpreter','none');
end

setfig('Total Ion Count');clf;
for i=1:length(testdata)
  tic(i)=sum(testdata{i}.TIC());
end
bar(tic);
set(gca,'XTick',1:length(testdata));
set(gca,'XTickLabel',cellfun(@(z) z.name, testdata,'Unif', false));
set(gca,'XTickLabelRotation',45);
set(gca,'ticklabelinterpreter','none');
ylabel('Total Ion Count');


return;

setfig('Aligned TIC');clf;
leg={};
for i=1:length(testdata)
  t=interp1(testmap{i}.time(:,2),testmap{i}.time(:,1),testdata{i}.time);
  testdata{i}.plotTIC(false,t);
  hold on;
  leg{end+1}=testdata{i}.name;
end
xlabel('Aligned time');
title('Aligned TIC');
l=legend(leg,'Location','best');
set(l,'Interpreter','none');
  

clist=[];
for i=1:length(testdata)
  testdata{i}.clusterpeaks('maxpeaks',50);
  c=testdata{i}.clusters;
  clist(i,round([c.mz]))=[c.elution];
end
for i=1:size(clist,2)
  if sum(clist(:,i)>0)>1
    fprintf('%d %s\n', i, sprintf('%4d ',clist(:,i)));
  end
end

  