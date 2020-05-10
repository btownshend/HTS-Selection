% Check composition of each M/S file
ref=mzdata{end};
for i=1:length(mzdata)-1
  compounds.plotComposition(mzdata{i},'ref',ref);
  pause(0.1)
end
