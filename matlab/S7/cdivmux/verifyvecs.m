% Verify various statistics on the given set of vectors
function verifyvecs(vecs,mass,mzminsep)
fprintf('Verify %d vectors over %d targets\n', size(vecs));
% Check coverage
fprintf('%d:%d vectors/target\n',min(sum(vecs,1)),max(sum(vecs,1)));
fprintf('%d:%d targets/vector\n',min(sum(vecs,2)),max(sum(vecs,2)));
if nargin<2
  return;
end
assert(size(vecs,2)==length(mass));
% Check individual separation
sep=false(size(vecs,2));
lastnsep=0;
for i=1:size(vecs,1)
  sep(vecs(i,:),~vecs(i,:))=true;
  sep(~vecs(i,:),vecs(i,:))=true;
  nsep=sum(sum(~sep)==1);
  if nsep>0 && (lastnsep==0 || nsep==size(vecs,2) || (nsep>lastnsep+50 && mod(i,20)==0))
    fprintf('Using %d vectors, %d targets are uniquely separable\n', i, nsep);
    lastnsep=nsep;
  end
  if nsep==size(vecs,2)
    break;
  end
end

% Check mass overlap
if nargin<3
  mzminsep=[0.003,0.006,0.010];
end
for sep=mzminsep
  d=nan(length(mass));   % minimum sep
  for i=1:length(mass)
    [~,dtmp]=aliased(mass(i),mass([1:i-1,i+1:end]),sep);
    d(i,[1:i-1,i+1:end])=dtmp;
  end
  nalias=sum(d<sep);
  fprintf('Targets are indistinguishable within M/Z of %.4f from an average of %.1f other target\n', sep, mean(nalias));
  
  nalias=[];
  for i=1:size(vecs,1)
    nalias(i)=sum(sum(d(vecs(i,:),vecs(i,:))<sep));
  end
  fprintf('Within vectors, an average of %.1f [%d-%d] targets are indistinguishable within M/Z of %.4f from another target\n', mean(nalias),min(nalias),max(nalias),sep);
end

end
