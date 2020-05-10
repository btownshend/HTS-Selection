% Assume SDF loaded
hits={'01A4','31C3','31D4','31E3','31E4','31F10','41D4','41D10',...
      '41E2','41E3','41E4','41F2','41F3','41F6','41F10','91C2',...
      '91C3','91C4','91D2','91D4','91E2','91E3','91E4','91F2',...
      '91F3','91F4','91H11','101D9','101D11','101E6','101E10','101F11',...
      '31C2','31C8','41E7','41F7'};
for i=1:length(hits)
  h=hits{i};
  if i==1
    sel=sdf.find(h);
  else
    sel=sel|sdf.find(h);
  end
end
fprintf('Selected %d hits\n', sum(sel));
setfig('Hits'); clf;
sdf.plot(find(sel));

