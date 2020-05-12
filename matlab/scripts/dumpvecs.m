% Dump layout into files for use by robot program
file=fopen('ops.csv','w');
fprintf(file,'srcplate\tdestplate\tsrcwell');
for i=1:nvecspertarget
  fprintf(file,'\tdestwell%d',i);
end
fprintf(file,'\n');
for i=1:ngroups
  for j=1:nplatespergroup
    for k=1:ntargetsperplate
      fprintf(file,'G%dP%d\tG%dV\t%s',i,j,i,wellname80(k));
      v=find(gvecs(i,:,(j-1)*ntargetsperplate+k));
      for m=1:length(v)
        fprintf(file,'\t%s',wellname96(v(m)));
      end
      fprintf(file,'\n');
    end
  end
end
fclose(file);


function s=wellname80(i)
  s=sprintf('%c%02d',char('A'+mod(i-1,8)),floor((i-1)/8)+2);
end

function s=wellname96(i)
  s=sprintf('%c%02d',char('A'+mod(i-1,8)),floor((i-1)/8)+1);
end

