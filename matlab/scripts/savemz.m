subdir='S7';
r=input([subdir,'?'],'s');
if r~='y'
  return;
end
mzdir=['/Users/bst/TMP/mzdata/',subdir];
fprintf('Saving mzdata in %s...',mzdir);
for i=1:length(mzdata)
  fprintf('%d...',i);
  mzd=mzdata{i};
  fname=sprintf('%s/%d.mat',mzdir,i);
  save(fname,'mzd');
end
fprintf('done\n');
clear mzd;
