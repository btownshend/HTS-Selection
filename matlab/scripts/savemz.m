dir='/Users/bst/TMP/mzdata';
fprintf('Saving mzdata in %s...',dir);
for i=1:length(mzdata)
  fprintf('%d...',i);
  mzd=mzdata{i};
  fname=sprintf('%s/%d.mat',dir,i);
  save(fname,'mzd');
end
fprintf('done\n');
clear mzd;
