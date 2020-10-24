% Compare old and new QTOF
% Assume mzdata{1} has the old QTOF for V256B-B3, mzdata{1} has the new
% s7compounds has the old analysis (with samples{48} = V256B-B3
samp=28;
% Find some good targets (has good normic, should be present, and we've isolated) for M+H
hit=s7compounds.normic(:,1,samp)>.5 & s7compounds.contains(:,samp) & isfinite(s7compounds.meantime);
% summary shows 129/130/256, 16 false positives for samp
fprintf('Have %d hits in %s\n', sum(hit), s7compounds.samples{samp});
fhit=find(hit);
mztol=[.003,.0005];
for ii=1:length(fhit)
  i=fhit(ii);
  mz=s7compounds.mass(i)+s7compounds.ADDUCTS(1).mass;
  s7compounds.getinfo(i);
  setfig(s7compounds.names{i});clf;
  subplot(211);
  mzdata{1}.ploteic(mz,'mztol',mztol(1),'newfig',false);
  ax=gca;
  subplot(212);
  mzdata{2}.ploteic(mz,'mztol',mztol(2),'newfig',false);
  ax(end+1)=gca;
  linkaxes(ax,'x');
  suptitle(s7compounds.names{i});
  pause
end

