% 101E6 seems to be a problem
ind=845;
compounds.getinfo(ind);
fprintf('\n');
truehits={'RowE','Col6','CDIV101','Full'};
falsehits={'RowD','RowG','RowH','Col2','Col5','Col9','CDIV061','CDIV081','CDIV111'};
all={truehits{:},falsehits{:}};
setfig('101E6');clf;
tiledlayout('flow');
mztarget=compounds.mztarget(ind);
mztol=.005;
for i=1:length(mzdata)
  if ismember(strrep(mzdata{i}.name,'.mzXML',''),all)
    [ic,mz,t]=mzdata{i}.mzscan(mztarget,'mztol',mztol);
    nexttile;
    plot(t,ic);
    ylabel('Ion Count');
    yyaxis right
    sel=ic>1e4;
    mzs=mz;
    mzs(~sel)=nan;   % Remove noisy measurements
    plot(t,mzs);
    ylabel('M/Z');
    ax=axis;
    ax(1:2)=[1400,1900];
    ax(3:4)=mztarget+mztol*[-1,1];
    axis(ax);
    title(sprintf('%s@%.4f',mzdata{i}.name,mztarget));
    fprintf('%s: m/z=%.4f+-%.4f\n', mzdata{i}.name, nanmean(mz), nanstd(mz));
  end
end
