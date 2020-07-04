% Generate and dump layout into files for use by robot program
wellsperplate=320;
v64plates={"V64A","V64B"};
vecsperplate=360;
targetspervec=64;
nvectors=vecsperplate*length(v64plates);
if ~exist('s7sdf','var')
  load('../../../data/matfiles/S7sdf.mat');
end
tmass=[s7sdf.sdf.MonoisotopicMass];
uplates=sort(unique({s7sdf.sdf.Plate384}));
topdown=true;

if topdown
  %v256=buildvecs(nvectors/4,tmass,targetspervec*4);
  v64parts={};
  vecs=false(nvectors,length(tmass));
  for i=1:size(v256,1)
    v64parts{i}=buildvecs(4,tmass(v256(i,:)),targetspervec);
    vecs((i-1)*4+(1:4),v256(i,:))=v64parts{i};
  end
else
  % Bottom-up
  vecs=buildvecs(vecsperplate*length(v64plates),tmass,targetspervec);
  v256=vecs(1:4:end,:)|vecs(2:4:end,:)|vecs(3:4:end,:)|vecs(4:4:end,:);
end

vecspertarget=unique(sum(vecs));
fprintf('Generated %d vectors from %d source plates with %d targets/vector, %d vectors/target\n', size(vecs,1), length(uplates), targetspervec, vecspertarget);

fprintf('V64:\n');
verifyvecs(vecs,tmass);

fprintf('\nV256:\n');
verifyvecs(v256,tmass);

file=fopen('V64.csv','w');
fprintf(file,'srcplate384\tsrcwell384\tsrcplate96\tsrcwell96\tmass');
for i=1:vecspertarget
  fprintf(file,'\tdestplate%d\tdestwell%d',i,i);
end
fprintf(file,'\n');
for i=1:length(s7sdf.sdf)
  s=s7sdf.sdf(i);
  fprintf(file,'%s\t%s\t%s\t%s\t%.4f',s.Plate384,s.Well384,s.BATCH_PLATE, s.BATCH_WELL, s.MonoisotopicMass);
  v=find(vecs(:,i));
  for m=1:length(v)
    fprintf(file,'\t%s\t%s',v64plates{ceil(v(m)/384)},wellname384(rem(v(m)-1,384)+1));
  end
  fprintf(file,'\n');
end
fclose(file);

function s=wellname384(i)
  s=sprintf('%c%02d',char('A'+mod(i-1,16)),floor((i-1)/16)+1);
end

