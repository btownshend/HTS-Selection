% Setup compounds with elution times for each compound
% Do this by distinguishing isomers by checking row, col, and plate totals separately -- as long as a compound is a unique isomer in at least one of these, then it should be distinguishable

datadir='../../data/';
matdir=[datadir,'matfiles/'];
resultsdir='../../results/';

if ~exist('s7sdf','var')
  load([matdir,'s7sdf.mat']);
end

if ~exist('s7vecs','var')
  s7vecs=load([matdir,'s7vecs.mat']);
end

mztol=.0005;

if ~exist('qsetup','var')
  qsetup=MSCompounds(mztol,40/60);
  qsetup.ADDUCTS=qsetup.ADDUCTS(1);   % Only M+H for now
  qsetup.addCompoundsFromSDF(s7sdf);
end


msdir2=[datadir,'MassSpec/20201002'];
mzmap2=[113,113;500,500];
timemap=[0,0;1,1];
data={ % Well, name, filename, mzmap, timemap,contains
       'A1','DMSO',[msdir2,'/A1.mzXML'],mzmap2,timemap,{},
       'B1','V5120',[msdir2,'/B1.mzXML'],mzmap2,timemap,{},
       'C1','V5120',[msdir2,'/C1.mzXML'],mzmap2,timemap,{},
       'D1','V2560A',[msdir2,'/D1.mzXML'],mzmap2,timemap,{},
       'E1','V2560B',[msdir2,'/E1.mzXML'],mzmap2,timemap,{},
       'F1','V2560A-sn',[msdir2,'/F1.mzXML'],mzmap2,timemap,{},
       'G1','V2560B-sn',[msdir2,'/G1.mzXML'],mzmap2,timemap,{},
       'H1','V2560A-pre',[msdir2,'/H1.mzXML'],mzmap2,timemap,{},
       'A2','V2560B-pre',[msdir2,'/A2.mzXML'],mzmap2,timemap,{},
       'E2','V256A-E6',[msdir2,'/E2.mzXML'],mzmap2,timemap,{},
       'F2','V256A-F6',[msdir2,'/F2.mzXML'],mzmap2,timemap,{},
       'G2','V256A-G6',[msdir2,'/G2.mzXML'],mzmap2,timemap,{},
       'H2','V256A-H6',[msdir2,'/H2.mzXML'],mzmap2,timemap,{},
          };
wells={'A3','B3','C3','D3','E3','F3','G3','H3',...
       'A4','B4','C4','D4','E4','F4','G4','H4',...
       'A5','B5','C5','D5','E5','F5','G5','H5',...
       'A6','B6','C6','D6'};
for i=1:length(wells)
  data(end+1,:)={wells{i},['V256A-',wells{i}],[msdir2,'/',wells{i},'.mzXML'],mzmap2,timemap,{}};
end
data=cell2struct(data,{'well','name','filename','mzmap','timemap','contains'},2);
for i=1:size(data,1)
  if isempty(data(i).contains)
    % Map from name to what the sample should contain
    data(i).contains=s7contains(qsetup,s7vecs,data(i).name);
  end
  if strcmp(data(i).name,'DMSO')
    data(i).conc=0;
  elseif strncmp(data(i).name,'V256A',5)
    data(i).conc=25e-9;
  elseif strcmp(data(i).well,'B1')
    data(i).conc=25e-9;
  else
    data(i).conc=100e-9;
  end
end

allmz=[];names={};
for i=1:length(qsetup.ADDUCTS)
  for j=1:length(qsetup.mass)
    allmz(end+1)=qsetup.mztarget(j,i);
    names{end+1}=[qsetup.names{j},'[',qsetup.ADDUCTS(i).name,']'];
  end
end

% Load data
if ~exist('mzdata','var')
  mzdata={};
end

for i=1:size(data,1)
  path=data(i).filename;
  if i>length(mzdata) || isempty(mzdata{i}) || ~strcmp(mzdata{i}.path,path)
      mztmp=MassSpec(path);
      mztmp.setLoad(data(i).conc*10e-6);
      % Prune out some data
      mztmp.filter([0,mztmp.time(end)-5],[min(allmz)-1,max(allmz)+1]);
                                                      %mztmp.keepmz(allmz,'mztol',mztol*2);
      mzdata{i}=mztmp;
  end
  mzdata{i}.name=data(i).name;   % Reset name to our list
  if ismember(mzdata{i}.path,qsetup.files)
    fprintf('Skipping reload of %s\n', mzdata{i}.name);
    continue;
  end

  % Build chromatograms if needed
  if isempty(mzdata{i}.featurelists) || ~any(strcmp({mzdata{i}.featurelists.src},'deconvolve'))
    if isempty(mzdata{i}.featurelists) || ~any(strcmp({mzdata{i}.featurelists.src},'buildchromatogram'))
      fprintf('Building chromatograms for %s...',mzdata{i}.name);
      mzdata{i}.targetedFeatureDetect(allmz, 'mztol',mztol,'names',names);
      fprintf('%d done\n',length(mzdata{i}.featurelists(end).features));
    end
    fprintf('Deconvolving chromatograms for %s...',mzdata{i}.name);
    mzdata{i}.deconvolve();
    fprintf('%d done\n',length(mzdata{i}.featurelists(end).features));
    %    fprintf('Finding isotopes for %s...',mzdata{i}.name);
    %    mzdata{i}.findisotopes();
    %    fprintf('done\n');
  end

  n=strsplit(data(i).name,'-');
  qsetup.addMS(mzdata{i},'group',n{1},'map',struct('mz',data(i).mzmap,'time',data(i).timemap),'contains',data(i).contains,'sample',data(i).name);
end

qsetup.findfeatures();

ref=3;
trace=[];
timetol=0.2;
usefiles=1:length(qsetup.samples);
usefiles=setdiff(usefiles,1:9);   % Skip pre/sn ones
minhits=2;
qsetup.assignTimes('clear',true,'timetol',timetol,'mztol',mztol,'trace',trace,'usefiles',usefiles,'minhits',minhits);
qsetup.checksensitivity(ref);
qsetup.assignTimes('clear',false,'timetol',timetol,'mztol',mztol,'trace',trace,'usefiles',usefiles,'minhits',minhits);
qsetup.assignTimes('clear',false,'timetol',timetol,'normicrange',[0.5,4.5],'trace',trace,'usefiles',usefiles,'minhits',minhits);
qsetup.assignTimes('clear',false,'timetol',timetol,'maxFP',1,'trace',trace,'usefiles',usefiles,'minhits',minhits);
qsetup.assignTimes('clear',false,'timetol',timetol,'maxFN',1,'trace',trace,'usefiles',usefiles,'minhits',minhits);
qsetup.assignTimes('clear',false,'timetol',timetol,'maxFN',1,'maxFP',1,'trace',trace,'usefiles',usefiles,'minhits',minhits);
qsetup.checksensitivity(ref);
qsetup.assignTimes('clear',false,'timetol',timetol,'maxFN',2,'maxFP',2,'trace',trace,'usefiles',usefiles,'minhits',minhits);


report=qsetup.report();

qsetup.checkmzoffset();
qsetup.checktime(ref,'timetol',qsetup.TIMEFUZZ/2);
qsetup.checksensitivity(ref);

qsetup.summary();

ds=datestr(now,'YYYYmmDDHHMM');
writetable(report,[resultsdir,sprintf('report-%s.csv',ds)]);
fprintf('Saving compounds...');
save(sprintf('%s/qsetup-%s.mat',matdir,ds),'qsetup');
fprintf('done\n');

for i=1:length(qsetup.ADDUCTS)
  qsetup.pcolorplot(i);
end

% Summarize by compound
namelist={'5A02'};
for i=1:length(namelist)
  qsetup.getinfo(namelist{i})
  fprintf('\n');
end
