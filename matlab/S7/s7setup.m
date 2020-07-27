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

if ~exist('compounds','var')
  compounds=Compounds();
  compounds.addCompoundsFromSDF(s7sdf);
end


msdir=[datadir,'MassSpec/'];
mzmap=[133,133;521,521];
timemap=[0,0;1,1];
data={ % Well, name, filename, mzmap, timemap,contains
          'A1','DMSO','?',mzmap,timemap,{},
          'B1','V256B-D3','?',mzmap,timemap,{},
          'C1','PR80-A1','?',mzmap,timemap,{},
          'D1','PR80-H11','?',mzmap,timemap,{},
          'E1','R320-A','?',mzmap,timemap,{},
          'F1','R320-P','?',mzmap,timemap,{},
          'G1','V5120','?',mzmap,timemap,{},
          'H1','CDIQ0005-Col3','?',mzmap,timemap,{'5A02','5B02','5C02','5D02','5E02','5F02','5G02','5H02','7A02','7B02','7C02','7D02','7E02','7F02','7G02','7H02'},
          'A2','CDIQ0005-Col22','?',mzmap,timemap,{'6A11','6B11','6C11','6D11','6E11','6F11','6G11','6H11','8A11','8B11','8C11','8D11','8E11','8F11','8G11','8H11'},
          'B2','DMSO','?',mzmap,timemap,{},
          'C2','V256A-A1','?',mzmap,timemap,{},
          'D2','V256A-B1','?',mzmap,timemap,{},
          'E2','V256A-C1','?',mzmap,timemap,{},
          'F2','V256A-D1','?',mzmap,timemap,{},
          'G2','V256A-E1','?',mzmap,timemap,{},
          'H2','V256A-F1','?',mzmap,timemap,{},
          'A3','V256A-G1','?',mzmap,timemap,{},
          'B3','V256A-H1','?',mzmap,timemap,{},
          'C3','V256A-A2','?',mzmap,timemap,{},
          'D3','V256A-B2','?',mzmap,timemap,{},
          'E3','V256A-C2','?',mzmap,timemap,{},
          'F3','V256A-D2','?',mzmap,timemap,{},
          'G3','V256A-E2','?',mzmap,timemap,{},
          'H3','V256A-F2','?',mzmap,timemap,{},
          'A4','V256A-G2','?',mzmap,timemap,{},
          'B4','V256A-H2','?',mzmap,timemap,{},
          'C4','V256A-A3','?',mzmap,timemap,{},
          'D4','V256A-B3','?',mzmap,timemap,{},
          'E4','V256A-C3','?',mzmap,timemap,{},
          'F4','V256A-D3','?',mzmap,timemap,{},
          'G4','V256B-A1','?',mzmap,timemap,{},
          'H4','V256B-B1','?',mzmap,timemap,{},
          'A5','V256B-C1','?',mzmap,timemap,{},
          'B5','V256B-D1','?',mzmap,timemap,{},
          'C5','V256B-E1','?',mzmap,timemap,{},
          'D5','V256B-F1','?',mzmap,timemap,{},
          'E5','V256B-G1','?',mzmap,timemap,{},
          'F5','V256B-H1','?',mzmap,timemap,{},
          'G5','V256B-A2','?',mzmap,timemap,{},
          'H5','V256B-B2','?',mzmap,timemap,{},
          'A6','V256B-C2','?',mzmap,timemap,{},
          'B6','V256B-D2','?',mzmap,timemap,{},
          'C6','V256B-E2','?',mzmap,timemap,{},
          'D6','V256B-F2','?',mzmap,timemap,{},
          'E6','V256B-G2','?',mzmap,timemap,{},
          'F6','V256B-H2','?',mzmap,timemap,{},
          'G6','V256B-A3','?',mzmap,timemap,{},
          'H6','V256B-B3','?',mzmap,timemap,{},
          'A7','V256B-C3','?',mzmap,timemap,{},
          'B7','DMSO','?',mzmap,timemap,{}
          };
data=cell2struct(data,{'well','name','filename','mzmap','timemap','contains'},2);
for i=1:size(data,1)
  if isempty(data(i).contains)
    % Map from name to what the sample should contain
    data(i).contains=s7contains(compounds,s7vecs,data(i).name);
  end
end

% Load data
if ~exist('mzdata','var')
  mzdata={};
end

if length(mzdata)>size(data,1)
  fprintf('Removing last %d entries from mzdata\n', length(mzdata)-size(data,1));
  mzdata=mzdata(1:size(data,1));
end

for i=1:size(data,1)
  path=[msdir,'/',data(i).filename];
  if i>length(mzdata) || ~strcmp(mzdata{i}.path,path)
      mzdata{i}=MassSpec(path);
      mzdata{i}.setLoad(2000*1e-15);
      % Prune out some data
      mzdata{i}.filter([0,mzdata{i}.time(end)-300],[130,530]);  % NOTE: this is using the localtimes to filter
  end
end

for i=1:size(data,1)
  n=strsplit(data(i).name,'-');
  s7compounds.addMS(mzdata{i},'group',n{1},'map',struct('mz',data(i).mzmap,'time',data(i).timemap),'contains',data(i).contains,'sample',data(i).name);
end

compounds.assignTimes();
compounds.summary();

report=compounds.report();

compounds.checkmzoffset();
ref=1;
compounds.checktime(ref,'timetol',compounds.TIMEFUZZ/2);
compounds.checksensitivity(ref);

writetable(report,[resultsdir,'report.csv']);
fprintf('Saving compounds...');
save([matdir,'compoundsS7.mat'],'compounds');
fprintf('done\n');

for i=1:length(compounds.ADDUCTS)
  compounds.pcolorplot(i);
end

% Summarize by compound
namelist={'5A02'};
for i=1:length(namelist)
  compounds.getinfo(namelist{i})
  fprintf('\n');
end
