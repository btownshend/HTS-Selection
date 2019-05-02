sdf=SDF();
% Load plates 1,11,...,111
sdf.read('ChemDivFull.sdf',111*80);
% Only keep the ones of interest
sdf=sdf.filter([1:10:111]);
sdf.getformulae();
sdf.csvwrite('dump.csv');
save('sdf.mat','sdf');
