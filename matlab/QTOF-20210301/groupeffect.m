ugroup=unique(v256.group);
tsens=nanmax(v256.tsens,[],2);
FN=arrayfun(@(z) nanmax([0,z.FN]),v256.astats);
hitgood=arrayfun(@(z) nanmax([0,z.hitgood]),v256.astats);
csel=FN==0;% & hitgood==3;
setfig('groupeffect');clf;
for i=1:length(ugroup)
  sel=find(strcmp(v256.group,ugroup{i}));
  sens=nan(length(v256.compound),length(sel));
  for jj=1:length(sel)
    j=sel(jj);
    sens(:,jj)=nanmax(v256.ic(:,:,j),[],2)./tsens;
    sens(~v256.contains(:,jj),jj)=nan;
  end
  nexttile;
  semilogy(v256.mass(csel),nanmean(sens(csel,:),2),'.');
  title(ugroup{i});
end