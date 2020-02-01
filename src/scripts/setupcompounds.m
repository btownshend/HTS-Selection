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


% Setup m/z, time remappings
maps=repmat(struct('mz',[],'time',[]),length(allfiles),1);
for i=1:length(allfiles)
  if ~isempty(strfind(allfiles(i).folder,'200124'))
    maps(i).mz=[133,133-10e-4;499,499-37e-4];
    maps(i).time=[0.3036    0.3476
                  2.0930    1.7762
                  2.5910    2.3867]*1e3;
  elseif ~isempty(strfind(allfiles(i).folder,'20190309'))
    maps(i).mz=[133,133-02e-4;499,499+21e-4];
    maps(i).time=[    0.3036    0.2811
                      2.5910    2.5782]*1e3;
  elseif ~isempty(strfind(allfiles(i).folder,'20190427'))
    maps(i).mz=[133,133-25e-4;499,499-62e-4];
    maps(i).time=[0 0; 1 1 ];
  elseif ~isempty(strfind(allfiles(i).folder,'20190501'))
    maps(i).mz=[133,133-4e-4;499,499+10e-4];
    maps(i).time=[0.3036    0.2835
                  2.4881    2.4791]*1e3;
  else
    error('No mapping for folder %s',allfiles(i).folder);
  end
end


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
      mzdata{i}=MassSpec(path);
      if ~isempty(strfind(allfiles{i},'200124'))
        mzdata{i}.setLoad(2000*1e-15);
      else
        mzdata{i}.setLoad(2500*1e-15);
      end
      % Prune out some data
      mzdata{i}.filter([300,2700],[127,505]);  % NOTE: this is using the localtimes to filter
  end
end


if ~exist('compounds','var')
  compounds=Compounds();
  compounds.addCompoundsFromSDF(sdf,'M+H');
  compounds.addCompoundsFromSDF(sdf,'M+Na');
end

for i=1:length(mzdata)
  if strncmp(mzdata{i}.name,'Full',4) || strncmp(mzdata{i}.name,'CDIV.',5) || strncmp(mzdata{i}.name,'CDIV-',5)
    compounds.addMS(mzdata{i},'group','Full','map',maps(i));
  elseif strncmp(mzdata{i}.name,'86',2)
    % (plates,columns) dropped out from 8619:8628 and 8630;  also the only contents of 8629
    dropouts=[1,4;31,3;31,4;31,7;31,8;31,9;41,6;41,7;51,8;61,4;91,2;91,3;91,4;91,5;101,6;101,8;101,10];
    assert(size(dropouts,1)==17);
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
      compounds.addMS(mzdata{i},'group','DiagPR','contains',contains,'map',maps(i));
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
      compounds.addMS(mzdata{i},'group','DiagPC','contains',contains,'map',maps(i));
    elseif id==8630
      for p=1:12
        for r=1:8
          for c=2:11
            if ~any((p-1)*10+1==dropouts(:,1) & c==dropouts(:,2))
              contains{end+1}=sprintf('%d%c%02d',(p-1)*10+1,r+'A'-1,c);
            end
          end
        end
      end
      assert(length(contains)==824);
      compounds.addMS(mzdata{i},'group','-Hits','map',maps(i),'contains',contains);
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
    compounds.addMS(mzdata{i},'contains',contains,'group','Col','map',maps(i));
  elseif strncmp(mzdata{i}.name,'Row',3)
    row=mzdata{i}.name(4);
    contains={};
    for p=1:12
      for c=2:11
        contains{end+1}=sprintf('%d%c%02d',(p-1)*10+1,row,c);
      end
    end
    compounds.addMS(mzdata{i},'contains',contains,'group','Row','map',maps(i));
  elseif strncmp(mzdata{i}.name,'CDIV',4)
    pnum=sscanf(mzdata{i}.name,'CDIV%d.mzXML');
    contains={};
    for r=1:8
      for c=2:11
        contains{end+1}=sprintf('%d%c%02d',pnum,r+'A'-1,c);
      end
    end
    compounds.addMS(mzdata{i},'contains',contains,'group','Plate','map',maps(i));
  elseif strncmp(mzdata{i}.name,'Full',4)
    compounds.addMS(mzdata{i},'group','Full','map',maps(i));
  elseif strcmp(mzdata{i}.name,'A2.mzXML')
    compounds.addMS(mzdata{i},'group','Individual','contains',{'31A2'},'map',maps(i));
  elseif strcmp(mzdata{i}.name,'A3.mzXML')
    compounds.addMS(mzdata{i},'group','Individual','contains',{'31A3'},'map',maps(i));
  else
    % Assume all compounds
    fprintf('Unable to decode filename "%s" -- ignoring\n',mzdata{i}.name);
  end
end

compounds.assignTimes();
compounds.summary();

report=compounds.report();

compounds.checkmzoffset();
ref=find(strcmp(compounds.files,'/Users/bst/Dropbox/SynBio/HTS-Selection/data/MassSpec/20190427 Row, Column/Full.mzXML'));
compounds.checktime(ref,compounds.TIMEFUZZ/2);

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
namelist={'1A04','91A02','91F03'};
for i=1:length(namelist)
  compounds.getinfo(namelist{i})
  fprintf('\n');
end
