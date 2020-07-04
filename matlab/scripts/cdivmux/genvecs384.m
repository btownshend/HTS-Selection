% Generate and dump layout into files for use by robot program
wellsperplate=320;
topdown=true;
v64plates={"V64A","V64B"};
v64perplate=360;
targetsperv64=64;
nv64=v64perplate*length(v64plates);

% Get target list from SDF
if ~exist('s7sdf','var')
  load('../../../data/matfiles/S7sdf.mat');
end
ntargets=length(s7sdf.sdf);
tmass=[s7sdf.sdf.MonoisotopicMass];

v64pertarget=nv64*targetsperv64/ntargets;
fprintf('Generating %d v64 vectors from %d source plates with %d targets/vector, %d vectors/target\n', nv64, length(unique({s7sdf.sdf.Plate384})), targetsperv64, v64pertarget);

if topdown
  % Build v256 vectors than split them for v64
  % Gives better v256 vectors
  fprintf('Build V256 vectors...\n');
  v256=buildvecs(nv64/4,tmass,targetsperv64*4);
  fprintf('Build V64 vectors...\n');
  v64parts={};
  v64=false(nv64,ntargets);
  for i=1:size(v256,1)
    v64parts{i}=buildvecs(4,tmass(v256(i,:)),targetsperv64);
    v64((i-1)*4+(1:4),v256(i,:))=v64parts{i};
  end
else
  % Bottom-up - build v64 vectors, then join them for v256
  fprintf('Build V64 vectors...\n');
  v64=buildvecs(v64perplate*length(v64plates),tmass,targetsperv64);
  fprintf('Build V256 vectors...\n');
  v256=v64(1:4:end,:)|v64(2:4:end,:)|v64(3:4:end,:)|v64(4:4:end,:);
end

fprintf('V64:\n');
verifyvecs(v64,tmass);

fprintf('\nV256:\n');
verifyvecs(v256,tmass);

file=fopen('V64.csv','w');
fprintf(file,'srcplate384\tsrcwell384\tsrcplate96\tsrcwell96\tmass');
for i=1:v64pertarget
  fprintf(file,'\tdestplate%d\tdestwell%d',i,i);
end
fprintf(file,'\n');
for i=1:length(s7sdf.sdf)
  s=s7sdf.sdf(i);
  fprintf(file,'%s\t%s\t%s\t%s\t%.4f',s.Plate384,s.Well384,s.BATCH_PLATE, s.BATCH_WELL, s.MonoisotopicMass);
  v=find(v64(:,i));
  for m=1:length(v)
    fprintf(file,'\t%s\t%s',v64plates{ceil(v(m)/384)},wellname384(rem(v(m)-1,384)+1));
  end
  fprintf(file,'\n');
end
fclose(file);
save('vecs.mat','v64','v256');

