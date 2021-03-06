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
%diagdata=[];
%for k=8632:8640
%  diagdata=[diagdata;dir([msdir,sprintf('200124-Diags-8630/%d.mzXML',k)])];
%end
diagdata=dir([msdir,'20200310/*.mzXML']);
diag2data=dir([msdir,'20200225/*.mzXML']);

allfiles=[rowdata;coldata;platedata;fulldata;indivdata;diag2data;diagdata];

%allfiles=allfiles(1:20);   % For testing

% Setup m/z, time remappings
maps=repmat(struct('mz',[],'time',[]),length(allfiles),1);
for i=1:length(allfiles)
  if ~isempty(strfind(allfiles(i).folder,'200124'))
    maps(i).mz=[133,133-7e-4;521,521-41e-4];
    maps(i).time=[304 346
                  2015 1700
                  2591 2380]/60;
  elseif ~isempty(strfind(allfiles(i).folder,'20200225'))
    maps(i).mz=[133,133+3e-4;521,521+20e-4];
    maps(i).time=[203 202
                  2894 2880]/60;
  elseif ~isempty(strfind(allfiles(i).folder,'20200310'))
    maps(i).mz=[133,133-9e-4;521,521-5e-4];
    maps(i).time=[203 211
                  2894 2891]/60;
  elseif ~isempty(strfind(allfiles(i).folder,'20190309'))
    maps(i).mz=[133,133-1e-4;521,521+20e-4];
    maps(i).time=[203 179
                  2894 2883]/60;
  elseif ~isempty(strfind(allfiles(i).folder,'20190427'))
    maps(i).mz=[133,133-24e-4;521,521-66e-4];
    maps(i).time=[0 0; 1 1 ];
  elseif ~isempty(strfind(allfiles(i).folder,'20190501'))
    maps(i).mz=[133,133+1e-4;521,521+6e-4];
    maps(i).time=[204 185
                  2894 2884]/60;
  else
    error('No mapping for folder %s',allfiles(i).folder);
  end
end


% Load data
if ~exist('mzdata','var')
  mzdata={};
end

if ~exist('compounds','var')
  compounds=MSCompounds();
  compounds.addCompoundsFromSDF(sdf);
  %  compounds.addCompoundsFromSDF(sdf.filter(sdf.find(1,[],2)));
end

for i=1:length(allfiles)
  % Can skip if already done (but always redo last file since it may have been incomplete)
  if i<length(compounds.files)
    fprintf('Skipping mzdata{%d} ... already done\n', i);
    continue;
  end

  % Load mzdata if needed
  path=[allfiles(i).folder,'/',allfiles(i).name];
  if i>length(mzdata) || isempty(mzdata{i}) || ~strcmp(mzdata{i}.path,path)
      mzdata{i}=MassSpec(path);
      if ~isempty(strfind(allfiles(i).folder,'200124'))
        mzdata{i}.setLoad(2000*1e-15);
      else
        mzdata{i}.setLoad(2500*1e-15);
      end
      % Prune out some data
      mzdata{i}.filter([200,2900]/60,[130,530]);  % NOTE: this is using the localtimes to filter
  end

  % Build chromatograms if needed
  if isempty(mzdata{i}.featurelists) || ~any(strcmp({mzdata{i}.featurelists.src},'deconvolve'))
    if isempty(mzdata{i}.featurelists) || ~any(strcmp({mzdata{i}.featurelists.src},'buildchromatogram'))
      fprintf('Building chromatograms for %s...',mzdata{i}.name);
      mzdata{i}.buildchromatograms();
      fprintf('%d done\n',length(mzdata{i}.featurelists(end).features));
    end
    fprintf('Deconvolving chromatograms for %s...',mzdata{i}.name);
    mzdata{i}.deconvolve();
    fprintf('%d done\n',length(mzdata{i}.featurelists(end).features));
    fprintf('Finding isotopes for %s...',mzdata{i}.name);
    mzdata{i}.findisotopes();
    fprintf('done\n');
  end
  
  if strncmp(mzdata{i}.name,'Full',4) || strncmp(mzdata{i}.name,'CDIV.',5) || strncmp(mzdata{i}.name,'CDIV-',5)  || strncmp(mzdata{i}.name,'old-CDIV',8)
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
  elseif strncmp(mzdata{i}.name,'DMSO',4)
    compounds.addMS(mzdata{i},'group','None','contains',{'NONE'},'map',maps(i));
  elseif strcmp(mzdata{i}.name,'101H6.mzXML')
    compounds.addMS(mzdata{i},'group','Individual','contains',{'101H06'},'map',maps(i));
  elseif strcmp(mzdata{i}.name,'31C2.mzXML')
    compounds.addMS(mzdata{i},'group','Individual','contains',{'31C02'},'map',maps(i));
  elseif strcmp(mzdata{i}.name,'41D7.mzXML')
    compounds.addMS(mzdata{i},'group','Individual','contains',{'41D07'},'map',maps(i));
  elseif strcmp(mzdata{i}.name,'91A3.mzXML')
    compounds.addMS(mzdata{i},'group','Individual','contains',{'91A03'},'map',maps(i));
  elseif strcmp(mzdata{i}.name,'91C2.mzXML')
    compounds.addMS(mzdata{i},'group','Individual','contains',{'91C02'},'map',maps(i));
  elseif strcmp(mzdata{i}.name,'A2.mzXML')
    compounds.addMS(mzdata{i},'group','Individual','contains',{'31A02'},'map',maps(i));
  elseif strcmp(mzdata{i}.name,'A3.mzXML')
    compounds.addMS(mzdata{i},'group','Individual','contains',{'31A03'},'map',maps(i));
  else
    % Assume all compounds
    fprintf('Unable to decode filename "%s" -- ignoring\n',mzdata{i}.name);
  end
end

compounds.findfeatures();
doassign;
compounds.checkmzoffset();
ref=find(strcmp(compounds.files,'/Users/bst/Dropbox/SynBio/HTS-Selection/data/MassSpec/20190427 Row, Column/Full.mzXML'));
compounds.checktime(ref,'timetol',compounds.TIMEFUZZ/2);
compounds.checksensitivity(ref);

compounds.summary();

report=compounds.report();
writetable(report,[resultsdir,'report.csv']);
%fprintf('Saving compounds...');
%save([matdir,'compounds.mat'],'compounds');
%fprintf('done\n');

if false
  fprintf('Saving mzdata...');
  save([matdir,'mzdata.mat'],'mzdata');
  fprintf('done\n');
end       

for i=1:length(compounds.ADDUCTS)
  compounds.pcolorplot(i);
end

% Summarize by compound
namelist={'1A04','91A02','91F03'};
for i=1:length(namelist)
  compounds.getinfo(namelist{i})
  fprintf('\n');
end
