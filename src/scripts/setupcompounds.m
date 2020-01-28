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
indivdata=[]; % dir([msdir,'20190506-Individual/*.mzXML']);
diagdata=dir([msdir,'200124-Diags-8630/*.mzXML']);

allfiles=[rowdata;coldata;platedata;fulldata;indivdata;diagdata];
mzoffsets=[repmat(.0028,1,length(rowdata)),repmat(.0025,1,length(coldata)),repmat(.0089,1,length(platedata)),.0038,repmat(0.0036,1,length(indivdata)),repmat(0.0043,1,length(diagdata))];
% Offsets of rerun on 5/1 seem different
mzoffsets(2)=.0076;
mzoffsets([10,18])=.0086;
% Setup time offsets
identityMap=[0 0; 1 1 ];

% Load data
if ~exist('mzdata','var')
  mzdata={};
end

if length(mzdata)>length(allfiles)
  fprintf('Removing last %d entries from mzdata\n', length(mzdata)-length(allfiles));
  mzdata=mzdata(1:length(allfiles));
end

for i=1:length(allfiles)
  path=[allfiles(i).folder,'/',allfiles(i).name];
  if i>length(mzdata) || ~strcmp(mzdata{i}.path,path)
      mzdata{i}=MassSpec(path,'mzoffset',mzoffsets(i));
      mzdata{i}.setLoad(2500*1e-15);
      % Prune out some data
      mzdata{i}.filter([300,2700],[127,505]);  % NOTE: this is using the localtimes to filter
  end
end

% Set time lookups
for i=1:length(mzdata)
  mzdata{i}.timeTLU=identityMap;
end
for i=length(mzdata)-length(diagdata)+1:length(mzdata)
  mzdata{i}.timeTLU=[0.3933    0.3669
                     1.7714    2.0839
                     2.5138    2.6997]*1e3;
end

% Set m/z offsets
for i=1:length(mzdata)
  mzdata{i}.adjmzoffset(mzoffsets(i));
end

if ~exist('compounds','var')
  compounds=Compounds();
  compounds.addCompoundsFromSDF(sdf,'H');
end

for reps=1:1
  fprintf('Pass %d...\n',reps);
  for i=1:length(mzdata)
    if strncmp(mzdata{i}.name,'Full',4) || strncmp(mzdata{i}.name,'CDIV.',5) || strncmp(mzdata{i}.name,'CDIV-',5)
      compounds.addMS(mzdata{i},'group','Full');
    elseif strncmp(mzdata{i}.name,'86',2)
      id=str2num(mzdata{i}.name(1:4));
      contains={};
      if id>=8611 && id<=8618
        % row-plate diagonals
        for p=1:12
          for r=1:8
            if mod(r-p,8)==(id-8611)
              for c=2:11
                contains{end+1}=sprintf('%d%c%02d',(p-1)*10+1,r+'A'-1,c);
              end
            end
          end
        end
        assert(length(contains)==120);
        compounds.addMS(mzdata{i},'group','DiagPR','contains',contains);
      elseif id>=8631 && id<=8640
        % col-plate diagonals
        contains={};
        for p=1:12
          for c=1:10
            if mod(c-p,10)==(id-8631)
              for r=1:8
                contains{end+1}=sprintf('%d%c%02d',(p-1)*10+1,r+'A'-1,c+1);
              end
            end
          end
        end
        assert(length(contains)==96);
        compounds.addMS(mzdata{i},'group','DiagPC','contains',contains);
      elseif id==8630
        % TODO: some components dropped out
        compounds.addMS(mzdata{i},'group','-Hits');
      else
        fprintf('Unable to decode filename "%s" -- ignoring\n',mzdata{i}.name);
      end
    elseif strncmp(mzdata{i}.name,'Col',3)
      cnum=sscanf(mzdata{i}.name,'Col%d.mzXML');
      contains={};
      for p=1:12
        for r=1:8
          contains{end+1}=sprintf('%d%c%02d',(p-1)*10+1,r+'A'-1,cnum);
        end
      end
      compounds.addMS(mzdata{i},'contains',contains,'group','Col');
    elseif strncmp(mzdata{i}.name,'Row',3)
      row=mzdata{i}.name(4);
      contains={};
      for p=1:12
        for c=2:11
          contains{end+1}=sprintf('%d%c%02d',(p-1)*10+1,row,c);
        end
      end
      compounds.addMS(mzdata{i},'contains',contains,'group','Row');
    elseif strncmp(mzdata{i}.name,'CDIV',4)
      pnum=sscanf(mzdata{i}.name,'CDIV%d.mzXML');
      contains={};
      for r=1:8
        for c=2:11
          contains{end+1}=sprintf('%d%c%02d',pnum,r+'A'-1,c);
        end
      end
      compounds.addMS(mzdata{i},'contains',contains,'group','Plate');
    elseif strncmp(mzdata{i}.name,'Full',4)
      compounds.addMS(mzdata{i},'group','Full');
    elseif strcmp(mzdata{i}.name,'A2.mzXML')
      compounds.addMS(mzdata{i},'group','Individual','contains',{'31A2'});
    elseif strcmp(mzdata{i}.name,'A3.mzXML')
      compounds.addMS(mzdata{i},'group','Individual','contains',{'31A3'});
    else
      % Assume all compounds
      fprintf('Unable to decode filename "%s" -- ignoring\n',mzdata{i}.name);
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
