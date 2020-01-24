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
fulldata=dir([msdir,'20190* Row, Column/Full.mzXML']);

allfiles=[rowdata;coldata;platedata;fulldata];
mzoffsets=[repmat(.0030,1,length(rowdata)),repmat(.0026,1,length(coldata)),repmat(.0087,1,length(platedata)),.0036];
% Offsets of rerun on 5/1 seem different
mzoffsets(2)=.0076;
mzoffsets([10,18])=.0086;
           
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

for i=1:length(mzdata)
  mzdata{i}.adjmzoffset(mzoffsets(i));
end

if ~exist('compounds','var')
  compounds=Compounds();
end


for reps=1:1
  fprintf('Pass %d...\n',reps);
  for i=1:length(mzdata)
    if strncmp(mzdata{i}.name,'Col',3)
      cnum=sscanf(mzdata{i}.name,'Col%d.mzXML');
      %      compounds.addFromSDF(mzdata{i},sdf.filter(sdf.find([],[],cnum)),'group','Col');
      compounds.addFromSDF(mzdata{i},sdf,'contains',sdf.find([],[],cnum),'group','Col');
    elseif strncmp(mzdata{i}.name,'Row',3)
      row=mzdata{i}.name(4);
      %compounds.addFromSDF(mzdata{i},sdf.filter(sdf.find([],row)),'group','Row');
      compounds.addFromSDF(mzdata{i},sdf,'contains',sdf.find([],row),'group','Row');
    elseif strncmp(mzdata{i}.name,'CDIV',4)
      pnum=sscanf(mzdata{i}.name,'CDIV%d.mzXML');
      %compounds.addFromSDF(mzdata{i},sdf.filter(sdf.find(pnum)),'group','Plate');
      compounds.addFromSDF(mzdata{i},sdf,'contains',sdf.find(pnum),'group','Plate');
    elseif strncmp(mzdata{i}.name,'Full',4)
      compounds.addFromSDF(mzdata{i},sdf,'group','Full');
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
fprintf('Saving compounds...');
save([matdir,'compounds.mat'],'compounds');
fprintf('done\n');

if false
  fprintf('Saving mzdata...');
  save([matdir,'mzdata.mat'],'mzdata');
  fprintf('done\n');
end       

% Summarize by compound
namelist={'1A4','91A2','91F3'};
for i=1:length(namelist)
  compounds.getinfo(namelist{i})
  fprintf('\n');
end
