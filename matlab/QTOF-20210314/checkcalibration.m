files={'20210314initial/A12-2.mzXML','H2',
       '20210314initial/A12-5.mzXML','H5',
       '20210314initial/A12-10.mzXML','H10',
       '20210314initial/A12-20.mzXML','H20',
       '20210314initial/B12-2.mzXML','L2',
       '20210314initial/B12-5.mzXML','L5',
       '20210314initial/blank.mzXML','Blank'};
if ~exist('mz','var')
  mz={};
  for i=1:length(files),
    mz{i}=MassSpec(['../../data/MassSpec/',files{i,1}]);
  end
end

if false
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
end
setfig('All Ref Masses');clf;
mz{7}.checkcalibration('refmasses',[121.050873,149.02332,322.048121,922.009798,1221.990637,1521.971475]);

% Load compounds inlcuded in V256-A1
masses=mysql('SELECT c.monoisotopicMass FROM compounds.compounds c, compounds.mixtures m, compounds.contents co WHERE c.compound=co.compound AND m.mixture=co.mixture AND m.name=''V256A-A1''');
mztgt=masses+1.007276;

peakic=[];
time=[];
for i=1:size(files,1)
  fprintf('Scanning %s...',files{i,2});
  for j=1:length(mztgt)
    [ic,m,t]=mz{i}.mzscan(mztgt(j),'mztol',.005,'timerange',[2,25]);
    [peakic(i,j),ppos]=max(ic);
    time(i,j)=t(ppos);
  end
  fprintf('done\n');
end
% Find best time for each
besttime=median(time);
% Rescan for correct time
timetol=0.2;
peaks1=peakic;
for i=1:size(files,1)
  fprintf('Rescanning %s...',files{i,2});
  for j=1:length(mztgt)
    if abs(time(i,j)-besttime(j))>timetol
      [ic,m,t]=mz{i}.mzscan(mztgt(j),'mztol',.005,'timerange',besttime(j)+timetol*[-1,1]);
      [peakic(i,j),ppos]=max(ic);
      time(i,j)=t(ppos);
    end
  end
  fprintf('done\n');
end

setfig('Peak IC');clf;
tl=tiledlayout('flow');
ref=1;
ax=[];
for i=setdiff(1:7,ref)
  nexttile;
  sel=abs(time(i,:)-time(ref,:))<timetol;  % Within 0.2 minute
  ratio=peakic(i,:)./peakic(ref,:);
  loglog(peakic(ref,sel),ratio(sel),'.');
  title(sprintf('%s %.2fx N=%d',files{i,2},nanmean(ratio(sel&peakic(ref,:)>1e4)),sum(sel)));
  ylabel(sprintf('%s/%s',files{i,2},files{ref,2}));
  ax(end+1)=gca;
end
xlabel(tl,files{ref,2});
title(tl,'m/z in V256A-A1 using [M+H]');
linkaxes(ax);
atmp=axis;
atmp(3:4)=[0.1/2,20];
axis(atmp);
