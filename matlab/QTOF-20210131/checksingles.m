% Check all the single compound mass spec runs
for i=1:length(mzdata)
  if ~strncmp(mzdata{i}.name,'CDIQ',4)
    continue;
  end
  mzd=mzdata{i};
  tgt=find(qsetup.contains(:,i));
  assert(length(tgt)==1);
  fprintf('Checking %s(%d) for %s (%d)\n', mzd.name, i, v256.names{tgt}, tgt);
  adduct=v256.astats(tgt).adduct;
  if isempty(adduct)
    fprintf('Not found in V256, plotting entire EIC for M+H\n');
    setfig(mzd.name);clf;
    tiledlayout('flow');
    for i=1:length(v256.ADDUCTS)
      nexttile;
      mzd.ploteic(v256.mztarget(tgt,i),'mztol',0.0005,'newfig',false);
    end
  else
    mzd.ploteic(v256.mztarget(tgt,adduct),'timerange',v256.timewindow(tgt,:)+[-.1,.1],'mztol',0.0005,'newfig',true);
    hold on;
    fprintf('Expected peak at m/z=%.5f, T=%.2f with intensity %.0f\n', v256.mztarget(tgt,adduct), v256.meantime(tgt),v256.tsens(tgt,adduct));
    plot(v256.meantime(tgt),v256.tsens(tgt,adduct),'*r','HandleVisibility','off');
    ax=axis;
    ax(4)=max(ax(4),v256.tsens(tgt,adduct)*1.1);
    axis(ax);
    plot(v256.timewindow(tgt,1)*[1,1],ax(3:4),':r','HandleVisibility','off');
    plot(v256.timewindow(tgt,2)*[1,1],ax(3:4),':r','HandleVisibility','off');
  end
end