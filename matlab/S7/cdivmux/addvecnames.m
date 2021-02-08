% When the s7vecs were built, the structure did not contain any indication of the ordering of the vectors
% add that to the struct
v64names={};
v64plates={"V64A","V64B"};
for i=1:720
  v64names{i}=sprintf('%s-%s',v64plates{ceil(i/384)},wellname384(rem(i-1,384)+1));
end
v256names={};
v256plates={"V256A","V256B"};
for i=1:180
  pos=mod(i-1,96)+1;
  col=floor((pos-1)/8)+1;
  row='A'+mod(pos-1,8);
  v256names{i}=sprintf('%s-%s%d',v256plates{ceil(i/96)},row,col);
end
s7vecs.v64names=v64names;
s7vecs.v256names=v256names;
s7vecs.cdiqnames=arrayfun(@(z) sprintf('%s-%s',s7vecs.csvdata.srcplate384{z},s7vecs.csvdata.srcwell384{z}),1:5120,'Unif',false);
s7vecs.cdivnames=arrayfun(@(z) sprintf('%s-%s',s7vecs.csvdata.srcplate96{z},s7vecs.csvdata.srcwell96{z}),1:5120,'Unif',false);
