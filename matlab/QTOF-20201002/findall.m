adducts=struct('name',{'M+H','M+Na','M+K'},'mass',{1.007276,22.989218,38.963158});

for i=10:length(mzdata)
  cont=find(qsetup.contains(:,i));
  for j=1:length(cont)
    mzdata{i}.findtarget(qsetup.sdf.sdf(cont(j)).Formula,'adducts',adducts,'mztol',0.002,'noise',500,'timetol',.6,'dbsave',true);
  end
end
