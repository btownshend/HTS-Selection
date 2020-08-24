n=length(allfl);
%fl=FeatureList.commonfeatures(allfl(1:n),'minintensity',2000,'minmatch',n-2,'timetol',1,'mztol',.015);
allmz=vertcat(fl.features.mzrange);
alltime=vertcat(fl.features.timerange);
mz=nanmedian(allmz')';
time=nanmean(alltime')';
for pnum=1:n
setfig('alignall');clf;
subplot(221);
mzrep=repmat(mz,1,size(allmz(:,pnum),2));
plot(mz,allmz(:,pnum)-mzrep,'.');
ppm=nanmedian((allmz(:,pnum)-mzrep)./mzrep*1e6);
fprintf('Delta m/z = %.1f PPM\n', ppm);
hold on;
ax=axis;
plot(ax(1:2),ax(1:2)*(ppm/1e6));
xlabel('m/z');
ylabel('Delta m/z');
title(allfl{pnum}.name);
subplot(222);
hist((allmz(:,pnum)-mzrep)./mzrep*1e6,100);

ax=axis;
ax(1)=min(ax(1),-100);
ax(2)=max(ax(2),100);
axis(ax);
xlabel('PPM change m/z');
ylabel('Count');
title(allfl{pnum}.name);
subplot(223);
plot(time,alltime(:,pnum)-repmat(time,1,size(alltime(:,pnum),2)),'.');
xlabel('time');
ylabel('Delta time');
input('continue?');
end