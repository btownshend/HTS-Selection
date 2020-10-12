% Setup compounds with elution times for each compound
% Do this by distinguishing isomers by checking row, col, and plate totals separately -- as long as a compound is a unique isomer in at least one of these, then it should be distinguishable

datadir='../../data/';
matdir=[datadir,'matfiles/'];
resultsdir='../../results/';


msdir=[datadir,'MassSpec/20201002'];
mzmap=[113,113;500,500];
timemap=[0,0;1,1];
data={ % Well, name, filename, mzmap, timemap,contains
          'A1','DMSO','A1.mzXML',mzmap,timemap,{},
          'B2','Acic200u','B2.mzXML',mzmap,timemap,{},
          'C2','Acic2m','C2.mzXML',mzmap,timemap,{},
          'D2','AcicFresh','D2.mzXML',mzmap,timemap,{},
          'E6','DMSO','E6.mzXML',mzmap,timemap,{},
       'B1','V5120','B1.mzXML',mzmap,timemap,{},
          };
data=cell2struct(data,{'well','name','filename','mzmap','timemap','contains'},2);

% Load data
if ~exist('mzdata','var')
  mzdata={};
end

acicmass=225.208;
for i=1:size(data,1)
  path=[msdir,'/',data(i).filename];
  if i>length(mzdata) || isempty(mzdata{i}) || ~strcmp(mzdata{i}.path,path)
      mzdata{i}=MassSpec(path);
      mzdata{i}.setLoad(100e-9*20e-6);
      % Prune out some data
      mzdata{i}.filter([0,mzdata{i}.time(end)-5],acicmass+[0,40]);  % NOTE: this is using the localtimes to filter
  end
  mzdata{i}.name=data(i).name;   % Reset name to our list

  continue;
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
end

