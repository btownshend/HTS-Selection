ds=datestr(now,'YYYYmmDDHHMM');
diary(['v256only-',ds,'.log']);
fsel=[];
for col=1:8
  for row='A':'H'
    if col==8 && row>='E'
      break;
    end
    vname=sprintf('V256A-%c%d',row,col);
    ind=find(strcmp(qsetup.samples,vname));
    if isempty(ind)
      fprintf('Missing vector: %s\n', vname);
    else
      fsel(end+1)=ind;
    end
  end
end
v256=qsetup.copypart('fsel',fsel);
v256.MZFUZZ=0.0012;
v256.TIMEFUZZ=0.2;

minhits=2;
trace=[];

v256.findfeatures();   % Should already exist from qsetup

% First run to be able to set fsens
v256.fsens(:)=1;
v256.assignTimes('clear',true,'normicrange',[0.1,10],'trace',trace,'minhits',minhits);
v256.checksensitivity();  % Initial fsens
% Start over with good fsens
v256.assignTimes('clear',true,'trace',trace,'minhits',minhits);
fprintf('Allow 1 false negative\n');
v256.assignTimes('clear',false,'trace',trace,'minhits',minhits,'maxFN',1);
fprintf('Allow 3 false positives\n');
v256.assignTimes('clear',false,'maxFP',3,'trace',trace,'minhits',minhits);
fprintf('Allow 3 false positives and 1 false negative\n');
v256.assignTimes('clear',false,'maxFP',3,'maxFN',1,'trace',trace,'minhits',minhits);
fprintf('Increased time tol\n');
v256.assignTimes('clear',false,'timetol',v256.TIMEFUZZ*2,'maxFP',0,'maxFN',0,'trace',trace,'minhits',minhits);
fprintf('Increased normicrange\n');
v256.assignTimes('clear',false,'normicrange',[0.2,5],'trace',trace,'minhits',minhits);
v256.checksensitivity(); % Final
v256.platesummary;

report=v256.report();
resultsdir='../../results/';
writetable(report,'/tmp/report.csv','QuoteStrings',true);
system(sprintf('sed -e "s/NaN//g" /tmp/report.csv > %s',[resultsdir,sprintf('report-%s.csv',ds)]));

fprintf('Saving compounds...');
allfeatures=v256.allfeatures;
reffeatures=v256.reffeatures;
v256.allfeatures=[];
v256.reffeatures=[];
save(sprintf('%s/v256-%s.mat',matdir,ds),'v256');
save(sprintf('%s/v256-allfeatures-%s.mat',matdir,ds),'allfeatures');
save(sprintf('%s/v256-reffeatures-%s.mat',matdir,ds),'reffeatures');
v256.allfeatures=allfeatures;
v256.reffeatures=reffeatures;

fprintf('done\n');
diary off


