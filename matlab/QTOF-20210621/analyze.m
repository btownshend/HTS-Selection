% Analyze ChemSpace comparisons
db=NGSDatabase.getinstance();
db.open();

if ~exist('chemcompare','var')
cmd=sprintf('SELECT r.msrun, s.name, r.name, r.mixture, r.well, r.concentration, r.injection FROM massspec.msruns r, massspec.mssessions s WHERE r.mssession=s.mssession AND s.mssession=6 AND r.well != ''''');
[msrun,sessionname, runname,mixture,well,conc,inject]=mysql(cmd);
msdata=struct('msrun',num2cell(msrun),'sessionname',sessionname,'runname',runname,'mixture',num2cell(mixture),'well',well,'conc',num2cell(conc),'inject',num2cell(inject));
well={msdata.well};
[~,ord]=sort(cellfun(@(z) [z(2:end),z(1)],well,'Unif',false));

chemcompare={};
for i=1:length(msdata)
  d=msdata(i);
  mix=Mixtures.instance().get(d.mixture);
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
  chemcompare{i}=mzd;
end
end

mztol=0.005;
adducts=[1.007276,22.989218,38.963158];
anames={'M+H','M+Na','M+K'};
linked=[];
prev=[];
for ia=1:length(adducts)
  setfig(sprintf('peaks - %s',anames{ia}));clf;
  tl=tiledlayout(7,2);
  for i=1:length(chemcompare)
    c=chemcompare{i};
    d=msdata(i);
    mix=Mixtures.instance().get(d.mixture);
    if isempty(mix.contents)
      continue;
    end
    co=Compounds.instance().get(mix.contents);
    fprintf('Well %s, %s, %s, %s, %.4f\n', d.well, d.runname, mix.name,co.name,co.mass)
    mz=co.mass+adducts(ia);
    nexttile;
    c.ploteic(mz,'mztol',mztol,'newfig',false);
    legend(co.name);
    if mz==prev
      linked(end+1)=gca;
      linkaxes(linked);
    else
      linked=gca;
      prev=mz;
    end
    
  end
  title(tl,anames{ia});
end

