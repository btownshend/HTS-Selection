% Remove a peak at a given M/Z from peak data and return ion count of peak
function [ppos,ic,ioncount,pkr]=removepeak(pks, time, mztarget, doplot)
  mztol=0.02;
  mzremovepeakwidth=[-0.1,3.1];
  timetol=60;  % Number of elution times to bracket peak (or if peak lost)
               % Find closest peak in each trace with mztol
  ioncount=zeros(length(pks),1);
  mzclosest=zeros(length(pks),1);
  for i=1:length(pks)
    [mzerr,mzclosest(i)]=min(abs(mztarget-pks{i}(:,1)));
    if mzerr<=mztol
      ioncount(i)=pks{i}(mzclosest(i),2);
    end
  end
  % Find elution trace with highest peak
  [maxic,tclosest]=max(ioncount);
  icthresh=maxic*.01;   % Cut off elutions at this fraction of max
  low=find(ioncount(1:tclosest-1)<icthresh,1,'last');
  low=max(low,tclosest-timetol);
  high=find(ioncount(tclosest+1:end)<icthresh,1)+tclosest;
  high=min(high,tclosest+timetol);
  ic=sum(ioncount(low:high));
  ppos=[pks{tclosest}(mzclosest(tclosest),1),low,tclosest,high];  % Position of peak [mz, tlow, tclosest, thigh]
  fprintf('Max peak for M/Z=%.4f found at M/Z=%.4f, T=[%4d-%4d-%4d], IC=%7.0f\n', mztarget, ppos, ic);
  % Remove this peak
  pkr=pks;
  for i=low:high
    sel=find((pks{i}(:,1)-mztarget)<=mzremovepeakwidth(2) & (pks{i}(:,1)-mztarget) >=mzremovepeakwidth(1));
    fprintf('Removed pks{%d}([%s])\n', i, sprintf('%d ',sel));
    pkr{i}(sel,2)=0;
  end
  % Plot if desired
  
  if nargin>=4 && doplot
    setfig(sprintf('removepeak %.4f',mztarget));clf;
    subplot(221);
    semilogy(time,ioncount);
    hold on;
    ax=axis;
    ax(1)=time(max(1,low-timetol));
    ax(2)=time(min(length(pks),high+timetol));
    axis(ax);
    plot(ax(1:2),icthresh*[1,1],':r');
    plot(time(tclosest)*[1,1],ax(3:4),':g');
    plot(time(low)*[1,1],ax(3:4),':r');
    plot(time(high)*[1,1],ax(3:4),':r');
    xlabel('Time');
    ylabel('Ion count');
    title(sprintf('M/Z %.4f (err=%.4f) IC=%.0f',mztarget,ppos(1)-mztarget,ic));

    subplot(223);
    stem(pks{tclosest}(:,1),pks{tclosest}(:,2));
    xlabel('m/z');
    ylabel('Ion count');
    title(sprintf('T(%d)=%f',tclosest,time(tclosest)));
  
    subplot(224);
    stem(pks{tclosest}(:,1),pkr{tclosest}(:,2));
    xlabel('m/z');
    ylabel('Ion count');
    title(sprintf('T(%d)=%f (after removal)',tclosest,time(tclosest)));
  end
end
