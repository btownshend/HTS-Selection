files={'20201002/blank.mzXML','blank',
       '20201002/A1.mzXML','DMSO',
       '20201002/A3.mzXML','V256A-A3',
       '2021.01.31/blank.mzXML','blank',
       '2021.01.31/A1.mzXML','DMSO',
       '2021.01.31/A3.mzXML','V256A-A2',
       '2021.03.01/blank1.mzXML','blank',
       '2021.03.01/A1.mzXML','DMSO',
       '2021.03.01/A3.mzXML','V256A-A2'};
if ~exist('mz','var')
  mz={};
  for i=1:length(files)
    mz{i}=MassSpec(['../../data/MassSpec/',files{i,1}]);
  end
end

setfig('Calibration');clf;
tl=tiledlayout('flow');
ax=[];
for i=1:length(files)
  nexttile;
  mz{i}.checkcalibration('newfig',false);
  ax(i)=gca;
  ti=get(gca,'title');
  ti.String=[ti.String,' ',files{i,2}];
  title(ti.String);
end
linkaxes(ax);



