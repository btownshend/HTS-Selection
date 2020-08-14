mzdir='/Users/bst/TMP/mzdata';
fprintf('Loading mzdata in %s...',mzdir);
for i=1:60
  fprintf('%d...',i);
  fname=sprintf('%s/%d.mat',mzdir,i);
  mzd=load(fname);
  mzdata{i}=mzd.mzd;
  mzdata{i}.featurelists=mzdata{i}.featurelists(1);
  mzdata{i}.mzxml=[];
end
fprintf('done\n');
clear mzd;