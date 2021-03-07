onlyc=[];
adducts=struct('name',{'M+H','M+Na','M+K','M+NH4','M+DMSO+H'},...
               'mass',{1.007276,22.989218,38.963158,18.033823,79.02122});  
% From https://fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator/
                                                                           %adducts=adducts(1);

debug=true;
diary findallhits.txt
fd=fopen('ident.csv','a');
fprintf(fd,'msrun,folder,file,compound,cname,adduct,m/z,time,intensity,nisotopes\n');
bestfeat=Feature.empty;bestmz=[];
cmd=sprintf('SELECT r.msrun, s.name, r.name, r.mixture, r.well, r.concentration, r.injection FROM massspec.msruns r, massspec.mssessions s WHERE r.mssession=s.mssession');
[msrun,sessionname, runname,mixture,well,conc,inject]=mysql(cmd);
compounds=Compounds();
compounds.load();
mixtures=Mixtures();
mixtures.load();
for i=1:length(msrun)
  mix=mixtures.get(mixture(i));
  if length(mix.contents)~=256
    continue;
  end
  if ~isempty(onlyc) && ~any(ismember(mix.contents,onlyc))
    continue;
  end
  tic
  fprintf('Processing msrun %s (%d) with mixture %s (%d)\n', name{i}, msrun(i), mix.name, mixture(i));
  path=sprintf('../../data/MassSpec/%s/%s.mzXML',sessionname{i},well{i});
  if ~exist(path,'file')
    error('File %s not found\n', path);
  end
  mzd=MassSpec(path);
  for j=1:length(mix.contents)
    c=compounds.get(mix.contents(j));
    if ~isempty(onlyc) && ~ismember(c.compound,onlyc)
      continue;
    end
    [bestfeature,bestmz,nisotopes]=mzd.findtarget(c.formula,'adducts',adducts,'mztol',0.002,'noise',500,'timetol',.6,'dbsave',false,'debug',debug);
    if isempty(bestfeature)
      fprintf(fd,'%d,"%s","%s",%d,"@%s",,,,,\n',msrun(i),sessionname{i},mzd.name,c.compound,c.name);
    else
      for k=1:length(bestfeature)
        fprintf(fd,'%d,"%s","%s",%d,"@%s"',msrun(i),sessionname{i},mzd.name,c.compound,c.name);
        fprintf(fd,',"%s",%.5f,%.2f,%.0f,%d\n',bestfeature(k).name,bestmz(k),bestfeature(k).time,bestfeature(k).intensity,nisotopes(k));
      end
    end
  end
  elapsed=toc;
  fprintf('Took %.0f seconds to process %s\n', name{i});
end
fclose(fd);
diary off
