minfold=3;
order=[];fold=[];
for i=2:length(summary.cleavage)
  c=summary.cleavage(i);
  if length(c.contains)==1 && c.conc==10e-6 && any(c.fold>minfold)
    order(end+1)=c.contains;
    fold(end+1)=max(c.fold);
  end
end
%order={'405D09','125F11','127E09','167A08','567C11','325H05','85B08','86B04','167H05','247E06'};
fd=fopen('order.csv','w');
fprintf(fd,'SMILES,name,fold\n');
for i=1:length(order)
  c=summary.compounds.get(order(i));
  fprintf(fd,'%s,"''%s",%.2f\n',c.smiles,c.name,fold(i));
end
fclose(fd);

    