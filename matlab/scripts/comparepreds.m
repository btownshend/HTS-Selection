% Assumes SDF is load
hits=readtable('~/Dropbox/SynBio/HTS-Selection/data/hits/Hits.csv');
pred=readtable('~/Dropbox/SynBio/HTS-Selection/python/moldist/src/preds.csv');
assert(length(hits.Properties.VariableNames)==length(pred.Properties.VariableNames));
assert(length(hits.Target)==length(pred.Target));
for i=2:length(hits.Properties.VariableNames)
  apt=hits.Properties.VariableNames{i};
  h=hits.(apt);
  p=pred.(apt);
  tp=h==1 & p==1;
  fp=h==0 & p==1;
  tn=h==0 & p==0;
  fn=h==1 & p==0;
  fprintf('%s: TP=%d, FP=%d, TN=%d, FN=%d\n', apt, sum(tp), sum(fp), sum(tn), sum(fn));

  setfig(apt);
  for i=1:3
    if i==1
      subplot(2,2,1);
      sdf.plot(tp);
      title('True Positives');
    elseif i==2
      subplot(2,2,2);
      sdf.plot(fn);
      title('False Negatives');
    else
      subplot(2,2,[3,4]);
      sdf.plot(fp,[],[],2);
      title('False Positives');
    end
  end
  suptitle(strrep(apt,'_','-'));
end
