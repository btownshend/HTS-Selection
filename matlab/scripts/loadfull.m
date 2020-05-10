% Load full (960 compound) runs

datadir='../../data/';
matdir=[datadir,'matfiles/'];
resultsdir='../../results/';

msdir=[datadir,'MassSpec/'];

files={'20190427 Row, Column/CDIV-full-conc-1.mzXML',
       '20190427 Row, Column/CDIV-full-conc-2.mzXML',
       '20190427 Row, Column/Full.mzXML',
       '20190327 Stability Experiments/CDIV_Lib_1.mzXML',
% '20190327 Stability Experiments/CDIV_Lib_2.mzXML',   % Bad pressure
       '20190327 Stability Experiments/CDIV_Lib_3.mzXML',
       '20190327 Stability Experiments/CDIV_Lib_4.mzXML',
       '190509-Round1TS/190509_T_T_T_8601.mzXML',
       '190509-Round1TS/190509_T_T_T_8602.mzXML',
       '190509-Round1TS/190509_T_T_T_8603.mzXML',
       '190509-Round1TS/190509_T_T_T_8604.mzXML',
       '190509-Round1TS/190509_T_T_T_8605.mzXML',
       '190509-Round1TS/190509_T_T_T_8606.mzXML',
       '20190327 Stability Experiments/Stab_25_1.mzXML',
       '20190327 Stability Experiments/Stab_25_2.mzXML',
       '20190327 Stability Experiments/Stab_25_3.mzXML',
       '20190327 Stability Experiments/Stab_25_4.mzXML',
       '20190327 Stability Experiments/Stab_30_1.mzXML',
       '20190327 Stability Experiments/Stab_30_2.mzXML',
       '20190327 Stability Experiments/Stab_30_3.mzXML',
%       '20190327 Stability Experiments/Stab_30_4.mzXML', % Bad Pressure
       '20190327 Stability Experiments/Stab_35_1.mzXML',
       '20190327 Stability Experiments/Stab_35_2.mzXML',
       '20190327 Stability Experiments/Stab_35_3.mzXML',
%       '20190327 Stability Experiments/Stab_35_4.mzXML', % Bad Pressure
       '20190327 Stability Experiments/Stab_40_1.mzXML',
%       '20190327 Stability Experiments/Stab_40_2.mzXML', % Bad Pressure
%       '20190327 Stability Experiments/Stab_40_3.mzXML', % Bad Pressure
       '20190327 Stability Experiments/Stab_40_4.mzXML',
       '20190309 CDiv Library Plates/CDIV Mix.mzXML'};

       
           
% Load data
if ~exist('mzfull','var')
  mzfull={};
end

mzoffset=0.0076;

for i=1:length(files)
  path=[msdir,'/',files{i}];
  if i>length(mzfull) || ~strcmp(mzfull{i}.path,path)
      mzfull{i}=MassSpec(path,'mzoffset',mzoffset);
      mzfull{i}.setLoad(2500*1e-15);
      % Prune out some data
      mzfull{i}.filter([300,2700],[127,505]);
  end
  fprintf('mzfull{%d} is %s from %s\n', i, mzfull{i}.name,mzfull{i}.path);
end

