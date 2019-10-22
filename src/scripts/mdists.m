% Setup molecular distances
datadir='../../data/';
matdir=[datadir,'matfiles/'];
resultsdir='../../results/';

if ~exist('sdf','var')
  load([matdir,'sdf.mat']);
end

m=MolDist();
m.setSDF(sdf);
m.loadcsv([datadir,'Tanimoto/dist20191021.csv']);
m.plotdist();
sets={'H960-751(720,072,594)',{'01A4','31E3','41F7','91D4','91E4','91F4','101E6'}, 
      'H960-506(251,228,172,561,875)',{'41E4','91C3','91D2','91E4','91F4'}, 
      'H960-650/003(319)',{'41F7','91D2','91E2','91F2'},
      'H960-050(724)',{'91C2','91C3','91C4'}, 
      'H969-113(892)',{'91C2','91D4','91E2','91E4','91F4'}, 
      'H960-940',{'91C3','91E4','91F2','91F4'}
     };
for i=1:size(sets,1)
  s=sets(i,:);
  x=m.closest(s{2},'mindist',0.0,'maxlist',20);
  m.plotclosest(x,s{1});
end
