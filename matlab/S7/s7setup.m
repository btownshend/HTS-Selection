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

if ~exist('s7compounds','var')
  s7compounds=Compounds();
  s7compounds.addCompoundsFromSDF(s7sdf);
end


msdir=[datadir,'MassSpec/20200724'];
mzmap=[113,113-17e-4;500,500-26e-4];
timemap=[0,0;1,1];
data={ % Well, name, filename, mzmap, timemap,contains
          'A1','DMSO','A1.mzXML',mzmap,timemap,{},
          'B1','V256B-D3','B1.mzXML',[113,113-27e-4;500,500-116e-4],timemap,{},
          'C1','PR80-A1','C1.mzXML',[113,113-26e-4;500,500-103e-4],timemap,{},
          'D1','PR80-H11','D1.mzXML',[113,113-20e-4;500,500-87e-4],timemap,{},
          'E1','R320-A','E1.mzXML',[113,113-20e-4;500,500-74e-4],timemap,{},
          'F1','R320-P','F1.mzXML',[113,113-21e-4;500,500-69e-4],timemap,{},
          'G1','V5120','G1.mzXML',[113,113-30e-4;500,500-45e-4],timemap,{},
          'H1','CDIQ0005-Col3','H1.mzXML',[113,113-22e-4;500,500-51e-4],timemap,{'5A02','5B02','5C02','5D02','5E02','5F02','5G02','5H02','7A02','7B02','7C02','7D02','7E02','7F02','7G02','7H02'},
          'A2','CDIQ0005-Col22','A2.mzXML',[113,113-18e-4;500,500-56e-4],timemap,{'6A11','6B11','6C11','6D11','6E11','6F11','6G11','6H11','8A11','8B11','8C11','8D11','8E11','8F11','8G11','8H11'},
          'B2','DMSO','B2.mzXML',mzmap,timemap,{},
          'C2','V256A-A1','C2.mzXML',mzmap,timemap,{},
          'D2','V256A-B1','D2.mzXML',mzmap,timemap,{},
          'E2','V256A-C1','E2.mzXML',mzmap,timemap,{},
          'F2','V256A-D1','F2.mzXML',mzmap,timemap,{},
          'G2','V256A-E1','G2.mzXML',mzmap,timemap,{},
          'H2','V256A-F1','H2.mzXML',mzmap,timemap,{},
          'A3','V256A-G1','A3.mzXML',mzmap,timemap,{},
          'B3','V256A-H1','B3.mzXML',mzmap,timemap,{},
          'C3','V256A-A2','C3.mzXML',mzmap,timemap,{},
          'D3','V256A-B2','D3.mzXML',mzmap,timemap,{},
          'E3','V256A-C2','E3.mzXML',mzmap,timemap,{},
          'F3','V256A-D2','F3.mzXML',mzmap,timemap,{},
          'G3','V256A-E2','G3.mzXML',mzmap,timemap,{},
          'H3','V256A-F2','H3.mzXML',mzmap,timemap,{},
          'A4','V256A-G2','A4.mzXML',mzmap,timemap,{},
          'B4','V256A-H2','B4.mzXML',mzmap,timemap,{},
          'C4','V256A-A3','C4.mzXML',mzmap,timemap,{},
          'D4','V256A-B3','D4.mzXML',mzmap,timemap,{},
          'E4','V256A-C3','E4.mzXML',mzmap,timemap,{},
          'F4','V256A-D3','F4.mzXML',mzmap,timemap,{},
          'G4','V256B-A1','G4.mzXML',mzmap,timemap,{},
          'H4','V256B-B1','H4.mzXML',mzmap,timemap,{},
          'A5','V256B-C1','A5.mzXML',mzmap,timemap,{},
          'B5','V256B-D1','B5.mzXML',mzmap,timemap,{},
          'C5','V256B-E1','C5.mzXML',mzmap,timemap,{},
          'D5','V256B-F1','D5.mzXML',mzmap,timemap,{},
          'E5','V256B-G1','E5.mzXML',mzmap,timemap,{},
          'F5','V256B-H1','F5.mzXML',mzmap,timemap,{},
          'G5','V256B-A2','G5.mzXML',mzmap,timemap,{},
          'H5','V256B-B2','H5.mzXML',mzmap,timemap,{},
          'A6','V256B-C2','A6.mzXML',mzmap,timemap,{},
          'B6','V256B-D2','B6.mzXML',mzmap,timemap,{},
          'C6','V256B-E2','C6.mzXML',mzmap,timemap,{},
          'D6','V256B-F2','D6.mzXML',mzmap,timemap,{},
          'E6','V256B-G2','E6.mzXML',mzmap,timemap,{},
          'F6','V256B-H2','F6.mzXML',mzmap,timemap,{},
          'G6','V256B-A3','G6.mzXML',mzmap,timemap,{},
          'H6','V256B-B3','H6.mzXML',mzmap,timemap,{},
% 'A7','V256B-C3','A7.mzXML',mzmap,timemap,{},
          'B7','DMSO','B7.mzXML',mzmap,timemap,{}
          };
data=cell2struct(data,{'well','name','filename','mzmap','timemap','contains'},2);
%data=data([3,5,6,7,8,2]);
for i=1:size(data,1)
  if isempty(data(i).contains)
    % Map from name to what the sample should contain
    data(i).contains=s7contains(s7compounds,s7vecs,data(i).name);
  end
end

% Load data
if ~exist('mzdata','var')
  mzdata={};
end

for i=1:size(data,1)
  path=[msdir,'/',data(i).filename];
  if i>length(mzdata) || isempty(mzdata{i}) || ~strcmp(mzdata{i}.path,path)
      mzdata{i}=MassSpec(path);
      mzdata{i}.setLoad(100e-9*20e-6);
      % Prune out some data
      mzdata{i}.filter([0,mzdata{i}.time(end)-5],[130,530]);  % NOTE: this is using the localtimes to filter
  end
  mzdata{i}.name=data(i).name;   % Reset name to our list
  if ismember(mzdata{i}.path,s7compounds.files)
    fprintf('Skipping reload of %s\n', mzdata{i}.name);
    continue;
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
    %    fprintf('Finding isotopes for %s...',mzdata{i}.name);
    %    mzdata{i}.findisotopes();
    %    fprintf('done\n');
  end

  n=strsplit(data(i).name,'-');
  s7compounds.addMS(mzdata{i},'group',n{1},'map',struct('mz',data(i).mzmap,'time',data(i).timemap),'contains',data(i).contains,'sample',data(i).name);
end

s7compounds.findfeatures();
doassign(s7compounds);

report=s7compounds.report();

s7compounds.checkmzoffset();
ref=7;
s7compounds.checktime(ref,'timetol',s7compounds.TIMEFUZZ/2);
s7compounds.checksensitivity(ref);

s7compounds.summary();

ds=datestr(now,'YYYYmmDDHHMM');
writetable(report,[resultsdir,sprintf('report-%s.csv',ds)]);
fprintf('Saving compounds...');
save(sprintf('%s/s7compounds-%s.mat',matdir,ds),'s7compounds');
fprintf('done\n');

for i=1:length(s7compounds.ADDUCTS)
  s7compounds.pcolorplot(i);
end

% Summarize by compound
namelist={'5A02'};
for i=1:length(namelist)
  s7compounds.getinfo(namelist{i})
  fprintf('\n');
end
