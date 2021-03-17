% Analyze CDIV565 rows
if ~exist('mz565','var')
cmd=sprintf('SELECT r.msrun, s.name, r.name, r.mixture, r.well, r.concentration, r.injection FROM massspec.msruns r, massspec.mssessions s WHERE r.mssession=s.mssession AND r.name LIKE ''CDIV0565-%%''');
[msrun,sessionname, runname,mixture,well,conc,inject]=mysql(cmd);
msdata=struct('msrun',num2cell(msrun),'sessionname',sessionname,'runname',runname,'mixture',num2cell(mixture),'well',well,'conc',num2cell(conc),'inject',num2cell(inject));

mz565={};
for i=1:length(msdata)
  d=msdata(i);
  fprintf('Processing msrun %s (%d) with mixture %s (%d)\n', d.runname, d.msrun, mix.name, d.mixture);
  path=sprintf('../../data/MassSpec/%s/%s.mzXML',d.sessionname,d.well);
  if ~exist(path,'file')
    error('File %s not found\n', path);
  end
  matfile=strrep(path,'.mzXML','.mat');
  if exist(matfile,'file')
    mzd=load(matfile);
    mzd=mzd.mztmp;
  else
    % Load mass spec from mzXML file
    mzd=MassSpec(path);
    mzd.moles=d.conc*d.inject;
    mztmp=mzd;
    save(matfile,'mztmp');
    clear mztmp;
  end
  mz565{i}=mzd;
end
end

mztol=0.005;
%un=MassSpec.uniquepeaks(mz565,'minic',2e5,'mztol',mztol);
cdiv=readtable('allmasses.csv');
allmass=cdiv.Var1;
allname=regexprep(cdiv.Var2,'CDIV0*(.*)-','$1');
allplates=str2double(regexprep(cdiv.Var2,'CDIV0*(.*)-.*','$1'));
pcnt=zeros(max(allplates),1);
contains=false(length(mz565),length(allmass));
for i=1:length(mz565)
  fprintf('%s: ',mz565{i}.name);
  mz=sort([un([un.maxind]==i).mz]);
  fprintf('%.4f ',mz);
  fprintf('\n');
  fprintf('M+H compounds: ');
  mass=mz-1.007276;
  for j=1:length(mass)
    ind=find(abs(mass(j)-allmass)<mztol);
    contains(i,ind)=true;
    pcnt(allplates(ind))=pcnt(allplates(ind))+1;
    ind564=ind(allplates(ind)==564);
    fprintf('%s ',strjoin(allname(ind564),','));
  end
  fprintf('\n');
end
[mp,ind]=max(pcnt);
fprintf('Plate %d has %d hits\n', ind, mp);
