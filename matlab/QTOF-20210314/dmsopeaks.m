% Analyze peaks in DMSO only by comparing with blank
%mzb=MassSpec('../../data/MassSpec/2021.03.14b/blank.mzXML');
%mzd=MassSpec('../../data/MassSpec/2021.03.14b/A1.mzXML');
%un=MassSpec.uniquepeaks({mzb,mzd},'minic',1e6);
ic=vertcat(un.ic);
% Find peaks that occur in DMSO but not in blank
sel=ic(:,2)>0 & ic(:,1)==0;
usel=un(sel);
[~,ord]=sort([usel.mz]);
usel=usel(ord);
% Group together isotopes
i=1;
while i<=length(usel)
  for j=i+1:length(usel)
    if abs(usel(j).mz-usel(i).mz(end)-1)<.01
      usel(i).mz=[usel(i).mz;usel(j).mz];
      usel(i).ic=[usel(i).ic;usel(j).ic];
      usel=usel([1:j-1,j+1:end]);
      i=i-1;
      break;
    end
  end
  i=i+1;
end

for i=1:length(usel)
  fprintf('%.4f ',usel(i).mz(1));
  for j=2:length(usel(i).mz)
    fprintf('+%.4f@%.2f ',usel(i).mz(j)-usel(i).mz(1),usel(i).ic(j,2)/usel(i).ic(1,2));
  end
  fprintf('\n');
end
