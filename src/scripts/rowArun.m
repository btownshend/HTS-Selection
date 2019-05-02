if ~exist('rowA','var')
  rowA=MassSpec('Row A P31.mzXML');
end
selsd=sdf.filter(31,'A');
searchPeaks=[selsd.sdf.MonoisotopicMass]+1+0.0068;

% Prune out some data
rowA.filter([400,2700],[min(searchPeaks)-5,max(searchPeaks)+5]);
%rowA.heatmap(searchPeaks);

% Find all the expected peaks
rowA.ident=[];
for i=1:length(searchPeaks)
  rowA.findcompound(searchPeaks(i),[],selsd.sdf(i));
end

for i=1:length(rowA.ident)
  rowA.plotident(i);
end

rowA.plotTIC();

return;

% Dot plot of peaks
%setfig('dotplot');clf;
msdotplot(peaks,time,'quantilevalue',0.9995);
hold on;
for i=1:length(searchPeaks)
  plot(searchPeaks(i),trange(end)+1,'rv');
end

if length(time)<100
  setfig('search');clf;
  plot3(repmat(mz,1,length(time)),repmat(time',n,1),y);
end

