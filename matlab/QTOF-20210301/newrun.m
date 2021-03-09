debug=true;

ds=datestr(now,'YYYYmmDDHHMM');
diary(['newrun-',ds,'.log']);

onlyc=3:5122;
adducts=struct('name',{'M+H','M+Na','M+K','M+NH4','M+DMSO+H'},...
               'mass',{1.007276,22.989218,38.963158,18.033823,79.02122});  
% From https://fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator/
                                                                           %adducts=adducts(1);

% Find all relevant mass spec runs
cmd=sprintf('SELECT r.msrun, s.name, r.name, r.mixture, r.well, r.concentration, r.injection FROM massspec.msruns r, massspec.mssessions s WHERE r.mssession=s.mssession');
[msrun,sessionname, runname,mixture,well,conc,inject]=mysql(cmd);
msdata=struct('msrun',num2cell(msrun),'sessionname',sessionname,'runname',runname,'mixture',num2cell(mixture),'well',well,'conc',num2cell(conc),'inject',num2cell(inject));
% Skip any from 2021.01.31 that were redone in 2021.03.01
mix3=mixture(strcmp(sessionname,'2021.03.01'));
skip=strcmp(sessionname,'2021.01.31') & ismember(mixture,mix3);
fprintf('Have %d/%d runs in 2021.01.31 that were later redone -- ignoring them\n', sum(skip), sum(strcmp(sessionname,'2021.01.31')));
msdata=msdata(~skip);
%assert(length(unique([msdata.mixture]))==length(msdata));   % Every mixture is unique
clear msrun sessionname runname mixture well conc inject

compounds=Compounds();
compounds.load();
mixtures=Mixtures();
mixtures.load();

v256=MSCompounds(.0005,40/60);
v256.MZFUZZ=0.0012;
v256.TIMEFUZZ=0.2;
v256.ADDUCTS=adducts;
falsenegs=onlyc([2:end,1]);
v256.addCompoundsFromDB(onlyc);

etotal=0;
for i=1:length(msdata)
  d=msdata(i);
  mix=mixtures.get(d.mixture);
  if length(mix.contents)~=256
    continue;
  end
  if ~isempty(onlyc) && ~any(ismember(mix.contents,onlyc))
    continue;
  end

  tic
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
  
  % Build a targetted feature list using expected masses plus an equal number of negatives (for balanced false positives)
  tgtmass=[];tgtnames={};
  for j=1:length(mix.contents)
    c=compounds.get(mix.contents(j));
    if ~isempty(onlyc) && ~ismember(c.compound,onlyc)
      continue;
    end
    tgtmass(end+1)=c.mass;
    tgtnames{end+1}=c.name;
    other=falsenegs(onlyc==mix.contents(j));
    c2=compounds.get(other);
    tgtmass(end+1)=c2.mass;
    tgtnames{end+1}=c2.name;
  end
  % Add all the adducts
  allmz=[];allnames={};
  for j=1:length(adducts)
    allmz=[allmz,tgtmass+adducts(j).mass];
    allnames=horzcat(allnames,cellfun(@(z) [z,'[',adducts(j).name,']'],tgtnames,'Unif',false));
  end
  
  fdmztol=2e-3;
  fprintf('Building chromatograms for %s against %d m/z with mztol=%f...',mzd.name,length(allmz),fdmztol);
  mzd.targetedFeatureDetect(allmz, 'mztol',fdmztol,'names',allnames,'debug',true);
  fprintf('%d done\n',length(mzd.featurelists(end).features));
  
  fprintf('Deconvolving chromatograms for %s...',mzd.name);
  mzd.deconvolve();
  fprintf('%d done\n',length(mzd.featurelists(end).features));

  v256.addMS(mzd,mix,'group',d.sessionname,'sample',d.runname);
  elapsed=toc;
  etotal=etotal+elapsed;
  fprintf('Processed %s in %.1f sec; %.0f sec to go\n', path, elapsed,etotal*(length(msdata)-i)/i);
end

v256.findfeatures();   

minhits=2;
maxfp=1;
trace=[];

% First run to be able to set fsens
v256.fsens(:)=1;
v256.assignTimes('clear',true,'normicrange',[0.1,10],'trace',trace,'minhits',minhits,'timetol',v256.TIMEFUZZ*4);
v256.checksensitivity();  % Initial fsens
% Start over with good fsens, but wide time range
v256.assignTimes('clear',true,'trace',trace,'minhits',minhits,'timetol',v256.TIMEFUZZ*4);
% Update times
v256.checktime('assign','dir','current',false);
v256.findfeatures();   
% Start over with good fsens, normal time range
v256.assignTimes('clear',true,'trace',trace,'minhits',minhits);
fprintf('Allow %d false positives\n',maxfp);
v256.assignTimes('clear',false,'maxFP',maxfp,'trace',trace,'minhits',minhits);
fprintf('Allow 1 false negative\n');
v256.assignTimes('clear',false,'trace',trace,'minhits',minhits,'maxFN',1);
fprintf('Allow %d false positives and 1 false negative\n',maxfp);
v256.assignTimes('clear',false,'maxFP',maxfp,'maxFN',1,'trace',trace,'minhits',minhits);
%fprintf('Increased time tol\n');
%v256.assignTimes('clear',false,'maxFP',maxfp,'maxFN',1,'timetol',v256.TIMEFUZZ*2,'trace',trace,'minhits',minhits);
fprintf('Increased mztol\n');
v256.assignTimes('clear',false,'maxFP',maxfp,'maxFN',1,'mztol',v256.MZFUZZ*2,'trace',trace,'minhits',minhits);
fprintf('Decreased mztol\n');
v256.assignTimes('clear',false,'maxFP',maxfp,'maxFN',1,'mztol',v256.MZFUZZ/2,'trace',trace,'minhits',minhits);
fprintf('Increased normicrange\n');
v256.assignTimes('clear',false,'maxFP',maxfp,'maxFN',1,'normicrange',[0.2,5],'trace',trace,'minhits',minhits);
v256.checksensitivity(); % Final

% Summarize run counts
rnum=arrayfun(@(z) max([z.run,0]),v256.astats);
setfig('AssignTimes');clf;
hist(rnum,0:max(rnum));
xlabel('AssignTimes Run');
ylabel('Compounds');
v256.platesummary;

report=v256.report();
resultsdir='../../results/';
for i=1:size(report)
  for j=1:length(report.file{i})
    report.(sprintf('file_%d',j)){i}=report.file{i}{j};
  end
end
report=removevars(report,'file');

writetable(report,'/tmp/report.csv','QuoteStrings',true);
system(sprintf('sed -e "s/NaN//g" /tmp/report.csv > %s',[resultsdir,sprintf('report-%s.csv',ds)]));

fprintf('Saving compounds...');
allfeatures=v256.allfeatures;
reffeatures=v256.reffeatures;
v256.allfeatures=[];
v256.reffeatures=[];
datadir='../../data/';
matdir=[datadir,'matfiles/'];
save(sprintf('%s/v256-%s.mat',matdir,ds),'v256');
save(sprintf('%s/v256-allfeatures-%s.mat',matdir,ds),'allfeatures');
save(sprintf('%s/v256-reffeatures-%s.mat',matdir,ds),'reffeatures');
v256.allfeatures=allfeatures;
v256.reffeatures=reffeatures;

fprintf('done\n');
diary off


