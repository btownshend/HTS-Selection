hts=HTS();
hts.load('SET_50000_pt1.csv');
hts.load('SET_50000_pt2.csv');
mw=nan(8,12,12);
for iplate=1:12;
  pname=sprintf('CDIV%04d',(iplate-1)*10+1);
  mw(:,:,iplate)=hts.getmw(pname);
end

% all plates
setfig('MW All');clf;
allmw=mw(:);allmw=allmw(isfinite(allmw));
fprintf('Have a total of %d MWs\n', length(allmw));
histogram(allmw,floor(min(allmw)):ceil(max(allmw)));
xlabel('MW');
ylabel('Count');
title(sprintf('All plates (N=%d)',length(allmw)));

% 31
setfig('MW 31');clf;
p31mw=mw(:,:,4);p31mw=p31mw(:);p31mw=p31mw(isfinite(p31mw));
fprintf('Have a total of %d MWs for plate 31\n', length(p31mw));
histogram(p31mw,floor(min(p31mw)):ceil(max(p31mw)));
xlabel('MW');
ylabel('Count');
title(sprintf('Plate CDIV0031 (N=%d)',length(p31mw)));

% 31 row A
setfig('MW 31-A');clf;
p31Amw=mw(1,:,4);p31Amw=p31Amw(:);p31Amw=p31Amw(isfinite(p31Amw));
fprintf('Have a total of %d MWs for plate 31, row A\n', length(p31Amw));
histogram(p31Amw,floor(min(p31Amw)):ceil(max(p31Amw)));
xlabel('MW');
ylabel('Count');
title(sprintf('Plate CDIV0031, Row A (N=%d)',length(p31Amw)));

