% Generate and dump layout into files for use by robot program
nplates=16;
wellsperplate=320;
v64plates={"V64A","V64B"};
vecsperplate=360;

%vecs=buildvecs(vecsperplate*length(v64plates),nplates*wellsperplate,64);
vecspertarget=unique(sum(vecs));
file=fopen('V64.csv','w');
fprintf(file,'srcplate\tsrcwell');
for i=1:vecspertarget
  fprintf(file,'\tdestplate%d\tdestwell%d',i,i);
end
fprintf(file,'\n');
for i=1:nplates
  for j=1:wellsperplate
    fprintf(file,'CDIV%d\t%s',i,wellname320(j));
    v=find(vecs(:,(i-1)*wellsperplate+j));
    for m=1:length(v)
      fprintf(file,'\t%s\t%s',v64plates{ceil(v(m)/384)},wellname384(rem(v(m)-1,384)+1));
    end
    fprintf(file,'\n');
  end
end
fclose(file);


function s=wellname320(i)
  s=sprintf('%c%02d',char('A'+mod(i-1,16)),floor((i-1)/16)+3);
end

function s=wellname384(i)
  s=sprintf('%c%02d',char('A'+mod(i-1,16)),floor((i-1)/16)+1);
end

