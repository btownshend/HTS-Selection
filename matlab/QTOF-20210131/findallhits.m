adducts=struct('name',{'M+H','M+Na','M+K','M+NH4','M+DMSO+H'},...
               'mass',{1.007276,22.989218,38.963158,18.033823,79.02122});  % From https://fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator/

debug=true;
!rm findallhits.txt 
diary findallhits.txt
fd=fopen('ident.csv','w');
fprintf(fd,'source,compound,cname,adduct,m/z,time,intensity,nisotopes\n');
for j=1:5120
  fprintf('Searching for compound %d: %s\n', j, qsetup.names{j});
  bestfeature={};
  for i=1:size(qsetup.contains,2)
    if sum(qsetup.contains(:,i))>256 || ~qsetup.contains(j,i)
      continue;
    end
    fprintf('Scanning %s (%d)\n',mzdata{i}.name,i);
    [bestfeature{i},bestmz,nisotopes]=mzdata{i}.findtarget(qsetup.sdf.sdf(j).Formula,'adducts',adducts,'mztol',0.002,'noise',500,'timetol',.6,'dbsave',false,'debug',debug);
    if ~isempty(bestmz)
      for k=1:length(bestmz)
        fprintf(fd,'%s,%d,"%s"',mzdata{i}.name,j,qsetup.names{j});
        fprintf(fd,',%s,%.5f,%.2f,%.0f,%d\n',bestfeature{i}(k).name,bestmz(k),bestfeature{i}(k).time,bestfeature{i}(k).intensity,nisotopes(k));
      end
    else
      fprintf(fd,'%s,%d,"%s"',mzdata{i}.name,j,qsetup.names{j});
      fprintf(fd,',,,,\n');
    end
  end
  % Find the best time for the entire group
  % TODO
end
diary off
