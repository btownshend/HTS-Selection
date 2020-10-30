adducts=struct('name',{'M+H','M+Na','M+K','M+NH4','M+DMSO+H'},...
               'mass',{1.007276,22.989218,38.963158,18.033823,79.02122});  % From https://fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator/

debug=true;
diary findallhits.txt
fd=fopen('ident.csv','a');
fprintf(fd,'source,compound,cname,adduct,m/z,time,intensity,nisotopes\n');
bestfeat=Feature.empty;bestmz=[];
for i=1:size(qsetup.contains,2)
  cont=find(qsetup.contains(:,i));
  if length(cont)>256
    continue;
  end
  for j=1:length(cont)
    [bestfeature,bestmz,nisotopes]=mzdata{i}.findtarget(qsetup.sdf.sdf(cont(j)).Formula,'adducts',adducts,'mztol',0.002,'noise',500,'timetol',.6,'dbsave',false,'debug',debug);
    fprintf(fd,'%s,%d,"%s"',mzdata{i}.name,cont(j),qsetup.names{cont(j)});
    if isfinite(bestmz)
      fprintf(fd,',%s,%.5f,%.2f,%.0f,%d\n',bestfeature.name,bestmz,bestfeature.time,bestfeature.intensity,nisotopes);
    else
      fprintf(fd,',,,,\n');
    end
  end
end
diary off
