% Find best M/Z peak in given time range
function mz=findmz(pks,time,range)
sel=find(time>=range(1) & time<=range(2));
if isempty(sel)
  error('No traces in range [%f,%f]\n', range);
end
highestic=0;
mz=[];
for ii=1:length(sel)
  i=sel(ii);
  [maxic,pos]=max(pks{i}(:,2));
  if maxic>highestic
    highestic=maxic;
    mz=pks{i}(pos,1);
  end
end
