x=readtable('ident.csv');
% Check if the same compound has consistent elution times
nmiss=0;nmultiple=0;nsingle=0;nbad=0;
etime=[];
for i=1:max(x.compound)
  xsel=x(x.compound==i,:);
  if range(xsel(isfinite(xsel.time),:).time)>1
    fprintf('Compound %d (%s) has inconsistent entries:\n',xsel.compound(1),xsel.cname{1});
    disp(xsel)
    nbad=nbad+1;
  elseif sum(isfinite(xsel.time))>1
    fprintf('Compound %d (%s) has %d consistent entries\n',xsel.compound(1),xsel.cname{1},height(xsel));
    nmultiple=nmultiple+1;
    etime(end+1)=nanmean(xsel.time);
  elseif any(isfinite(xsel.time))
    nsingle=nsingle+1;
  else
    nmiss=nmiss+1;
  end
end
fprintf('Have %d single-occurence entries, %d multiple consistent ones, %d inconsistent ones, %d missing\n', nsingle, nmultiple, nbad, nmiss);

% Plot elution time distribution using multiples only
setfig('Elution Time');clf;
histogram(etime,50);
xlabel('Elution Time (min)');
title('Elution time distribution for multiples');
