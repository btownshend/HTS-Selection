% Assume setup already run, so mzdata and qsetup are filled in
allmz=[];names={};
for i=1:length(qsetup.ADDUCTS)
  for j=1:length(qsetup.mass)
    allmz(end+1)=qsetup.mztarget(j,i);
    names{end+1}=[qsetup.names{j},'[',qsetup.ADDUCTS(i).name,']'];
  end
end
fprintf('Have %d target m/z\n', length(allmz));

% Use an mztol here based on resolving power; that gives us all the peaks that might include the desired peaks
% Assume an average mass of 300 and resolving power of 20000
fdmztol=300/2000;

% Build chromatograms if needed
for i=1:length(mzdata)
  force=length(mzdata{i}.featurelists)>=2 && mzdata{i}.featurelists(end-1).params.mztol~=fdmztol;
  if force|| isempty(mzdata{i}.featurelists) || ~any(strcmp({mzdata{i}.featurelists.src},'deconvolve'))
    tic
    if force || isempty(mzdata{i}.featurelists) || ~any(strcmp({mzdata{i}.featurelists.src},'buildchromatogram'))
      fprintf('Building chromatograms for %s against %d m/z with mztol=%f...',mzdata{i}.name,length(allmz),fdmztol);
      mzdata{i}.targetedFeatureDetect(allmz, 'mztol',fdmztol,'names',names);
      fprintf('%d done\n',length(mzdata{i}.featurelists(end).features));
    end
    fprintf('Deconvolving chromatograms for %s...',mzdata{i}.name);
    mzdata{i}.deconvolve();
    fprintf('%d done\n',length(mzdata{i}.featurelists(end).features));
    %    fprintf('Finding isotopes for %s...',mzdata{i}.name);
    %    mzdata{i}.findisotopes();
    %    fprintf('done\n');
    matpath=strrep(mzdata{i}.path,'.mzXML','.mat');
    mztmp=mzdata{i};
    fprintf('Saving matfile in %s...', matpath);
    save(matpath,'mztmp');
    fprintf('done\n');
    clear mztmp;
    toc
  end

  fl=mzdata{i}.featurelists(end);
  qsetup.allfeatures(i)=fl;
end

qsetup.findfeatures();

ref=3;
trace=[];
timetol=0.2;
usefiles=1:length(qsetup.samples);
usefiles=setdiff(usefiles,1:9);   % Skip pre/sn ones
minhits=2;
qsetup.assignTimes('clear',true,'timetol',timetol,'mztol',mztol,'trace',trace,'usefiles',usefiles,'minhits',minhits);
qsetup.checksensitivity(ref);
qsetup.assignTimes('clear',false,'timetol',timetol,'mztol',mztol,'trace',trace,'usefiles',usefiles,'minhits',minhits);
qsetup.assignTimes('clear',false,'timetol',timetol,'normicrange',[0.5,4.5],'trace',trace,'usefiles',usefiles,'minhits',minhits);
qsetup.assignTimes('clear',false,'timetol',timetol,'maxFP',1,'trace',trace,'usefiles',usefiles,'minhits',minhits);
qsetup.assignTimes('clear',false,'timetol',timetol,'maxFN',1,'trace',trace,'usefiles',usefiles,'minhits',minhits);
qsetup.assignTimes('clear',false,'timetol',timetol,'maxFN',1,'maxFP',1,'trace',trace,'usefiles',usefiles,'minhits',minhits);
qsetup.checksensitivity(ref);
qsetup.assignTimes('clear',false,'timetol',timetol,'maxFN',2,'maxFP',2,'trace',trace,'usefiles',usefiles,'minhits',minhits);


report=qsetup.report();

qsetup.checkmzoffset();
qsetup.checktime(ref,'timetol',qsetup.TIMEFUZZ/2);
qsetup.checksensitivity(ref);

qsetup.summary();

ds=datestr(now,'YYYYmmDDHHMM');
resultsdir='../../results/';
writetable(report,[resultsdir,sprintf('report-%s.csv',ds)]);
fprintf('Saving compounds...');
save(sprintf('%s/qsetup-%s.mat',matdir,ds),'qsetup');
fprintf('done\n');

for i=1:length(qsetup.ADDUCTS)
  qsetup.pcolorplot(i);
end

% Summarize by compound
namelist={'5A02'};
for i=1:length(namelist)
  qsetup.getinfo(namelist{i})
  fprintf('\n');
end
