if ~exist('blank','var')
  blank=MassSpec('Blank.mzXML');
end
% Prune out some data
%blank.filter([400,2700],[min(searchPeaks)-5,max(searchPeaks)+5]);
%blank.heatmap(searchPeaks);

dmso=78.0139355;
searchPeaks=[dmso];

% Find all the expected peaks
for i=1:length(searchPeaks)
  rowA.findcompound(searchPeaks(i));
  rowA.plotident(i);
end
blank.plotTIC();

