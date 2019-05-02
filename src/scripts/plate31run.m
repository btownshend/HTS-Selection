selsd=sdf.filter(31);
searchPeaks=[selsd.sdf.MonoisotopicMass]+1+0.0068;
searchNames={selsd.sdf.BATCH_WELL};

if ~exist('plate31','var')
  plate31=MassSpec('P31.mzXML');
  % Prune out some data
  plate31.filter([300,2700],[min(searchPeaks)-5,max(searchPeaks)+5]);
end
plate31.heatmap(searchPeaks);

% Find all the expected peaks
plate31.ident=[];
for i=1:length(searchPeaks)
  plate31.findcompound(searchPeaks(i),[], selsd.sdf(i));
end

for i=1:min(5,length(plate31.ident))
  plate31.plotident(i);
end

%color code for plate
pc=reshape([plate31.ident.ictotal],10,8);
pc(end+1,:)=nan;
pc(:,end+1)=nan;
setfig('Hits');clf;
pcolor(log10(pc)');
set(gca,'XTick',1.5:10.5);
set(gca,'XTickLabel',arrayfun(@(z) sprintf('%d',z),2:11,'UniformOutput',false));
set(gca,'YTick',1.5:8.5);
set(gca,'YTickLabel',{'A','B','C','D','E','F','G','H'});
axis ij;
colorbar;

plate31.plotTIC();

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

