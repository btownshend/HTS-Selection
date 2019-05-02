if ~exist('sdf','var')
  load ../HTBCFiles/sdf.mat
end

QQQBase='../QQQ/';

if ~exist('rowA','var')
  platedir=[QQQBase,'20190215 QTOF Plate 31/'];
  rowA=fullsetrun([platedir,'/Row A P31.mzXML'],8000);
  p31=fullsetrun([platedir,'/P31.mzXML'],1000);
end

if ~exist('stab1','var')
  stabdir=[QQQBase,'20190306 Stability Tests'];
  d=dir([stabdir,'/*.mzXML']);
  stab1={};
  for i=1:length(d)
    if isempty(strfind(d(i).name,'DMSO'))
      fmoles=0;
    else
      fmoles=833;
    end
    stab1{i}=fullsetrun(sprintf('%s/%s',stabdir,d(i).name),fmoles);
  end
end

if ~exist('plates','var')
  platedir=[QQQBase,'20190309 CDiv Library Plates'];
  plates={};
  for i=1:12
    plates{i}=fullsetrun(sprintf('%s/CDIV%03d.mzXML',platedir,i*10-9),2500);
  end
  mixed=fullsetrun(sprintf('%s/CDIV Mix.mzXML',platedir),417);
end

if ~exist('rc','var')
  rcdir=[QQQBase,'20190427 Row, Column'];
  d=dir([rcdir,'/*.mzXML']);
  rc={};
  for i=1:length(d)
    if isempty(strfind(d(i).name,'DMSO'))
      fmoles=0;
    else
      fmoles=833;
    end
    rc{i}=fullsetrun(sprintf('%s/%s',rcdir,d(i).name),fmoles);
  end
end

if ~exist('stab','var')
  stabdir=[QQQBase,'20190327 Stability Experiments'];
  d=dir([stabdir,'/*.mzXML']);
  stab={};
  for i=1:length(d)
    if isempty(strfind(d(i).name,'DMSO'))
      fmoles=0;
    else
      fmoles=208;
    end
    stab{i}=fullsetrun(sprintf('%s/%s',stabdir,d(i).name),fmoles);
  end
end

compounds=Compounds();
mzoffset=7.4e-3;   % Seems like all our QTOF M/Z are off by this much
compounds.addFromSDF(rowA,sdf.filter(31,'A'),'mzoffset',mzoffset);
for i=1:length(rc)
  if strncmp(rc{i}.name,'Col',3)
    cnum=sscanf(rc{i}.name,'Col%d.mzXML');
    compounds.addFromSDF(rc{i},sdf.filter([],[],cnum),'mzoffset',mzoffset);
  elseif strncmp(rc{i}.name,'Row',3)
    row=rc{i}.name(4);
    compounds.addFromSDF(rc{i},sdf.filter([],row),'mzoffset',mzoffset);
  elseif strncmp(rc{i}.name,'DMSO',4)
    ;
  else
    % Assume all compounds
    compounds.addFromSDF(rc{i},sdf,'mzoffset',mzoffset);
  end
end
compounds.addFromSDF(p31,sdf.filter(31),'mzoffset',mzoffset);
for i=1:length(plates)
  compounds.addFromSDF(plates{i},sdf.filter(i*10-9),'mzoffset',mzoffset);
end
compounds.addFromSDF(mixed,sdf,'mzoffset',mzoffset);
for i=1:length(stab1)
  compounds.addFromSDF(stab1{i},sdf,'mzoffset',mzoffset);
end
for i=1:length(stab)
  compounds.addFromSDF(stab{i},sdf,'mzoffset',mzoffset);
end
compounds.summary();
report=compounds.report();
fprintf('M/Z offset = %.4f\n', mzoffset+nanmedian(report.mzoffset));
writetable(report,'report.csv');

% Compute relative scaling
ref=15;
fprintf('Scaling relative to %s (%d)\n', compounds.files{ref},ref);
s=[];moles=[];
for i=1:length(compounds.files)
  s(i)=compounds.getscaling(ref,i);
  moles(i)=compounds.moles(i);
  fprintf('%2d %20.20s %5.0f %5.2f\n', i, compounds.files{i}, moles(i), s(i));
end
setfig('moles vs ic');clf;
loglog(moles*1e15,s,'o');
fit=sum(s)/sum(moles*1e15);
ax=axis;
hold on;
plot(ax(1:2),ax(1:2)*fit,':');
xlabel('Amount loaded (fmoles)');
ylabel('Relative Ion Count');
title('Ion Counts vs amount loaded');


function [obj,id]=fullsetrun(file,fmoles)
  obj=MassSpec(file);
  obj.setLoad(fmoles*1e-15);
  % Prune out some data
  obj.filter([300,2700],[127,505]);
  % Find all the expected peaks
  %obj.clusterpeaks('maxpeaks',5000);
end
