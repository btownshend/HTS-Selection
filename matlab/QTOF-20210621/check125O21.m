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
