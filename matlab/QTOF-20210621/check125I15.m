% Check results for 125I15 vs 8746
outliers=[565352879,565476652,565491775,566218011];   % naseqs that have no switching for 8746, but >2x for CDIQ
% x=summary108.plotgradsummary(ngs108,'metric','ratio_of_ratio','minswitching',1.1,'mintargets',1,'sort','tree','onlyconc',[]);
allsel2=[];
for i=1:length(outliers)
  sel=x.naseq==outliers(i);
  sel2=find(x.switching(:,sel)>=2);
  fprintf('Naseq %d has >=2-fold change for %d conditions\n', outliers(i), length(sel2));
  for j=1:length(sel2)
    ename=x.exptnames{sel2(j)};
    if ename(1)>='1' && ename(1)<='9'
      fprintf('%20.20s: %5.2f\n', ename, x.switching(sel2(j),sel))
      allsel2=union(allsel2,sel2(j));
    end
  end
end

% Check summary data
e1=find(strcmp({summary108.cleavage.target},'8746'));
e2=find(strcmp({summary108.cleavage.target},'CDIQ0125-I15'));
f1=summary108.cleavage(e1).fold;
f2=summary108.cleavage(e2).fold;

for j=1:length(allsel2)
  % Check other aptamers sensitive to this target
  sel3=find(x.switching(allsel2(j),:)>=2);
  sw=x.switching(allsel2(j),sel3);
  [~,ord]=sort(sw,'desc');
  sel3=sel3(ord);
  fprintf('%-20.20s:\n',x.exptnames{allsel2(j)});
  fprintf('      NASeq  Fold (8746/125I15)\n');
  for k=1:length(sel3)
    if true % ~ismember(x.naseq(sel3(k)),outliers)
      sel4=find(summary108.naseq==x.naseq(sel3(k)));
      fprintf('  %d %5.2f (%5.2f/%5.2f)\n',x.naseq(sel3(k)), x.switching(allsel2(j),sel3(k)),f1(sel4),f2(sel4))
    end
  end
end


return

% Check if CDIQ125O21 (CDIV125H11) may have been contaminated with CDIV125F11
% 125F11: M+H 425.0244 12.93min ic=19699
c=[]
c(end+1)=find(strcmp({msdata.runname},'8749'));
c(end+1)=find(strcmp({msdata.runname},'125O21'));
c(end+1)=find(strcmp({msdata.runname},'125K21'));
setfig('125O21 check');clf;
tiledlayout('flow');
ax=[];
for i=1:length(c)
  nexttile;
  chemcompare{c(i)}.ploteic(425.1333,'newfig',false);
  title(msdata(c(i)).runname);
  ax(end+1)=gca;
end
linkaxes(ax);
