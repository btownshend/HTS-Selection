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

if ~exist('qtc','var')
  qtc=Compounds();
  qtc.addCompoundsFromSDF(s7sdf);
end


msdir1=[datadir,'MassSpec/20200724'];
msdir2=[datadir,'MassSpec/20201002'];
mzmap1=[113,113-17e-4;500,500-26e-4];
mzmap2=[113,113;500,500];
timemap=[0,0;1,1];
data={ % Well, name, filename, mzmap, timemap,contains
       %'1.A1','DMSO',[msdir1,'/A1.mzXML'],mzmap1,timemap,{},
       %  '2.A1','DMSO',[msdir2,'/A1.mzXML'],mzmap2,timemap,{},
       %   '1.G6','V256B-A3',[msdir1,'/G6.mzXML'],mzmap1,timemap,{},
       '1.D4','V256A-B3',[msdir1,'/D4.mzXML'],mzmap1,timemap,{},
%  '1.H6','V256B-B3',[msdir1,'/H6.mzXML'],mzmap1,timemap,{},
%          '2.A3','V256B-A3',[msdir2,'/A3.mzXML'],mzmap2,timemap,{},
         '2.B3','V256A-B3',[msdir2,'/B3.mzXML'],mzmap2,timemap,{}
       
          };
data=cell2struct(data,{'well','name','filename','mzmap','timemap','contains'},2);
for i=1:size(data,1)
  if isempty(data(i).contains)
    % Map from name to what the sample should contain
    data(i).contains=s7contains(qtc,s7vecs,data(i).name);
  end
end

allmz=[];names={};
for i=1:1 % length(qtc.ADDUCTS)
  for j=1:length(qtc.mass)
    if s7compounds.contains(j,28)
      allmz(end+1)=qtc.mztarget(j,i);
      names{end+1}=[qtc.names{j},'[',qtc.ADDUCTS(i).name,']'];
    end
  end
end

% Load data
if ~exist('mzdata','var')
  mzdata={};
end

for i=1:size(data,1)
  if data(i).well(1) == '1'
    mztol=.01;
  else
    mztol=.0005;
  end
  path=data(i).filename;
  if i>length(mzdata) || isempty(mzdata{i}) || ~strcmp(mzdata{i}.path,path)
      mztmp=MassSpec(path);
      mztmp.setLoad(100e-9*20e-6);
      % Prune out some data
      mztmp.filter([0,mztmp.time(end)-5],[130,530]);  % NOTE: this is using the localtimes to filter
                                                      %mztmp.keepmz(allmz,'mztol',mztol*2);
      mzdata{i}=mztmp;
  end
  mzdata{i}.name=data(i).name;   % Reset name to our list
  if ismember(mzdata{i}.path,qtc.files)
    fprintf('Skipping reload of %s\n', mzdata{i}.name);
    continue;
  end

  % Build chromatograms if needed
  if isempty(mzdata{i}.featurelists) || ~any(strcmp({mzdata{i}.featurelists.src},'deconvolve'))
    if isempty(mzdata{i}.featurelists) || ~any(strcmp({mzdata{i}.featurelists.src},'targetedFeatureDetect'))
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
  qtc.addMS(mzdata{i},'group',n{1},'map',struct('mz',data(i).mzmap,'time',data(i).timemap),'contains',data(i).contains,'sample',data(i).name);
end

qtc.findfeatures();
return;

doassign;

report=qtc.report();

qtc.checkmzoffset();
ref=7;
qtc.checktime(ref,'timetol',qtc.TIMEFUZZ/2);
qtc.checksensitivity(ref);

qtc.summary();

ds=datestr(now,'YYYYmmDDHHMM');
writetable(report,[resultsdir,sprintf('report-%s.csv',ds)]);
fprintf('Saving compounds...');
save(sprintf('%s/qtc-%s.mat',matdir,ds),'qtc');
fprintf('done\n');

for i=1:length(qtc.ADDUCTS)
  qtc.pcolorplot(i);
end

% Summarize by compound
namelist={'5A02'};
for i=1:length(namelist)
  qtc.getinfo(namelist{i})
  fprintf('\n');
end
