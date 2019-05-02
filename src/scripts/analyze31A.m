if ~exist('mzxml','var')
  mzxml=mzxmlread('Row A P31.mzXML');
end
res=0.01;
rowApeaks=[306.1117,382.1133,424.1859,430.1311,286.0397,294.1256,214.0630,340.1099,327.1504,354.1692]+1;
%rowApeaks=rowApeaks(5);
[peaks,time]=mzxml2peaks(mzxml);
% Restrict to given time range
%trange=time([1,end]);
trange=[400,2700];
sel=find(time>=trange(1) & time<=trange(2));
time=time(sel);
peaks=peaks(sel);

mzrange=[floor(min(rowApeaks)-2),ceil(max(rowApeaks)+2)];
n=round(diff(mzrange)/res+1);
[mz,y]=msppresample(peaks,n,'Range',mzrange,'FWHH',res);
setfig('heatmap');clf;
msheatmap(mz,time,log(y),'resolution',res,'markers',rowApeaks);


% Find closest peaks
besttrace=[];
setfig('best by M/Z');clf;
setfig('best by T');clf;
nc=ceil(sqrt(length(rowApeaks)));
nr=ceil(length(rowApeaks)/nc);
for i=1:length(rowApeaks)
  height=[];
  [~,closest]=min(abs(mz-rowApeaks(i)));
  [level,besttrace(i)]=max(sum(y(closest+(-2:2),:),1));
  % Find actual peak corresponding to this
  [~,pp]=min(abs(peaks{besttrace(i)}(:,1)-rowApeaks(i)));
  mzmatch(i)=peaks{besttrace(i)}(pp,1);

  setfig('best by M/Z');
  subplot(nr,nc,i);
  plot(mz,y(:,besttrace(i)));
  xlabel('M/Z');
  hold on;
  ax=axis;
  plot(rowApeaks(i)*[1,1],ax(3:4),':r');
  title(sprintf('M/Z=%.4f(%.4f), T=%.2f, I=%.0f',rowApeaks(i),mzmatch(i),time(besttrace(i)),level));

  setfig('best by T');
  subplot(nr,nc,i);
  plot(time,y(closest,:));
  hold on;
  ax=axis;
  plot(time(besttrace(i))*[1,1],ax(3:4),':r');
  title(sprintf('M/Z=%.4f(%.4f), T=%.2f, I=%.0f',rowApeaks(i),mzmatch(i),time(besttrace(i)),level));
end

fprintf('Mean m/z error = %.4f +- %.4f\n', mean(mzmatch-rowApeaks), std(mzmatch-rowApeaks));

% Use new method of removing peaks
pr=peaks;
for i=1:length(rowApeaks)
  [ppos,ic,allic,pr]=removepeak(pr,time,rowApeaks(i),1);
end
% Resample without the identified peaks
[mzr,yr]=msppresample(pr,n,'Range',mzrange,'FWHH',res);

% Plot TIC showing best traces
setfig('TIC');clf;
plot(time,sum(y));
hold on;
plot(time,sum(yr));
title('Total Ion Chromatogram (TIC)')
xlabel('Retention Time')
ylabel('Relative Intensity')
legend('Original','Post Peak Removal');
hold on;
ax=axis;
for i=1:length(besttrace)
  plot(time(besttrace(i))*[1,1],[0,ax(4)*1.1],':r')
end

% Dot plot of peaks
%setfig('dotplot');clf;
msdotplot(peaks,time,'quantilevalue',0.9995);
hold on;
for i=1:length(rowApeaks)
  plot(rowApeaks(i),trange(end)+1,'rv');
end

if length(time)<100
  setfig('search');clf;
  plot3(repmat(mz,1,length(time)),repmat(time',n,1),y);
end

