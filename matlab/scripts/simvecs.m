% Simulate a situation with nhits and measurement using nvecs
% Return set of possible hits
function [pos,certain]=simvecs(vecs,nhits,nvectors)
  ntargets=size(vecs,2);
  if nargin<3
    nvectors=size(vecs,1);
  else
    vecs=vecs(1:nvectors,:);
  end
  
  hits=randsample(ntargets,nhits);
  hot=any(vecs(:,hits),2);
  pos=~any(vecs(~hot,:),1);
  assert(all(pos(hits)));
  certain=any(vecs(sum(vecs(:,pos),2)==1,:),1) & pos;   % Ones that have a unique target within the possible
  assert(all(ismember(find(certain),hits)));
  if sum(certain)>1 && sum(certain)<sum(pos)
    %keyboard;
    ;
  end
end
