if ~exist('itbl','var')
  il={};
  for i=1:length(mzdata)
    il{i}=mzdata{i}.featurelists(end).findisotopes();
    il{i}(:,end+1)=i;
  end
  ilist=vertcat(il{:});
  itbl=array2table(ilist,'VariableNames',{'f1','f2','mz','dmz','dt','ratio','ord1','ord2','file'});
end
setfig('Isotope statistics');clf;
tiledlayout('flow');

nexttile; % setfig('isotope m/z');clf;
hist(itbl.dmz-1.0030,200);
xlabel('\Delta m/z - 1.0030');

nexttile; % setfig('isotope time');clf;
hist(itbl.dt,100);
xlabel('\Delta T');

nexttile; % setfig('isotope ratio');clf;
leg={};
for ord=1:10
  sel=itbl.ord1==ord;
  if ~any(sel)
    break;
  end
  plot(itbl.mz(sel), itbl.ratio(sel)./itbl.mz(sel),'.');
  hold on;
  leg{end+1}=sprintf('Order %d',ord);
end
legend(leg);
set(gca,'YScale','log');
xlabel('m/z');
ylabel('Intensity ratio/(m/z)');

nexttile; % setfig('isotope order');clf;
boxplot(itbl.ratio./itbl.mz,itbl.ord1);
set(gca,'YScale','log');
xlabel('Order');
ylabel('Intensity ratio/(m/z)');
