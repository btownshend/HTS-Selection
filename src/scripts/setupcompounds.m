% Setup compounds with elution times for each compound
% Do this by distinguishing isomers by checking row, col, and plate totals separately -- as long as a compound is a unique isomer in at least one of these, then it should be distinguishable

datadir='../../data/';
matdir=[datadir,'matfiles/'];
resultsdir='../../results/';

if ~exist('sdf','var')
  load([matdir,'sdf.mat']);
end

msdir=[datadir,'MassSpec/'];

rowdata=dir([msdir,'20190* Row, Column/Row*.mzXML']);
rowdata=rowdata([1,9,3:8]);   % Use rerun
coldata=dir([msdir,'20190* Row, Column/Col*.mzXML']);
coldata=coldata([2,11,4:9,1,10]);   % Use rerun, reorder
platedata=dir([msdir,'20190309 CDiv Library Plates/CDIV*1.mzXML']);

allfiles=[rowdata;coldata;platedata];
mzoffsets=[repmat(-.0027,1,length(rowdata)),repmat(-.0027,1,length(coldata)),repmat(.0014,1,length(platedata))]+7.4e-3;
% Load data
if ~exist('mzdata','var')
  mzdata={};
end

for i=1:length(allfiles)
  path=[allfiles(i).folder,'/',allfiles(i).name];
  if i>length(mzdata) || ~strcmp(mzdata{i}.path,path)
      mzdata{i}=MassSpec(path,'mzoffset',mzoffsets(i));
      mzdata{i}.setLoad(2500*1e-15);
      % Prune out some data
      mzdata{i}.filter([300,2700],[127,505]);
  end
end


if ~exist('compounds','var')
  compounds=Compounds();
end


for reps=1:2
  fprintf('Pass %d...',reps);
  for i=1:length(mzdata)
    if strncmp(mzdata{i}.name,'Col',3)
      cnum=sscanf(mzdata{i}.name,'Col%d.mzXML');
      compounds.addFromSDF(mzdata{i},sdf.filter([],[],cnum),'group','Col');
    elseif strncmp(mzdata{i}.name,'Row',3)
      row=mzdata{i}.name(4);
      compounds.addFromSDF(mzdata{i},sdf.filter([],row),'group','Row');
    elseif strncmp(mzdata{i}.name,'CDIV',4)
      pnum=sscanf(mzdata{i}.name,'CDIV%d.mzXML');
      compounds.addFromSDF(mzdata{i},sdf.filter(pnum),'group','Plate');
    else
      % Assume all compounds
      fprintf('Unable to decode filename "%s" -- ignoring\n',mzdata{i});
    end
  end
end
compounds.summary();

report=compounds.report();

compounds.checkmzoffset();

writetable(report,[resultsdir,'report.csv']);
save([matdir,'compounds.mat'],'compounds');

% TODO, summarize by compound

