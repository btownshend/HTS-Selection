% Build a logic matrix size (nvectors,ntargets) such that each row has targetspervec true values
% Total number of trues in each column should be equal within 1
function [vecs,vecspertarget]=buildvecs(nvectors,mass,targetspervec)
  debug=1;
  ntargets=length(mass);
  vecs=false(nvectors,ntargets);
  vecspertarget=targetspervec*nvectors/ntargets;
  vecsperset=nvectors/vecspertarget;
    
  % Fill in set by set
  vecs=false(0,ntargets);
  for i=1:vecspertarget
    fprintf('%d ',i);

    if i>1
      % Compute matrix of number of distinguishing vectors between each pair of targets with vectors so far
      ndiffs=zeros(size(vecs,2));
      for i=1:size(vecs,1)
        ndiffs(vecs(i,:),~vecs(i,:))=ndiffs(vecs(i,:),~vecs(i,:))+1;
        ndiffs(~vecs(i,:),vecs(i,:))=ndiffs(~vecs(i,:),vecs(i,:))+1;
      end
      mindiffs=prctile(ndiffs(:),10);
      fprintf('Over vecs 1..%d, each target is distinguishable from an average of %.0f others in at least %.0f different vectors\n', ...
              size(vecs,1), mean(sum(ndiffs>=mindiffs)), mindiffs);
      vecs=[vecs;buildset(mass,vecsperset,debug,ndiffs)];
    else
      vecs=[vecs;buildset(mass,vecsperset,debug)];
    end
    if debug
      fprintf('\n');
    end
    
  end
  assert(all(sum(vecs)<=ceil(vecspertarget) & sum(vecs)>=floor(vecspertarget)));
  assert(all(sum(vecs,2)==targetspervec));
end

% Build a set of vectors where each target occurs exactly once and no two targets in the same vector have aliased masses
% Also, attempt to put targets that are indistinguishable based on 'indist' matrix into separate vectors, if possible (just a hint)
function vecs=buildset(mass,nvectors,debug,ndiffs)
  if nargin<4
    ndiffs=[];
  end
  targetspervec=length(mass)/nvectors;
  assert(floor(targetspervec)==targetspervec);   % Perfect fit
  vecs=false(nvectors,length(mass));
  mzminsep=[0.01,0.008,0.006,0.005]; 
  ndebug=0;
  cnt=zeros(length(mass),2,length(mzminsep));
  while sum(vecs(:))<length(mass)
    % Choose a target
    tbest=find(sum(vecs,1)==0);
    t=tbest(randi(length(tbest),1,1));
    valias=false(nvectors,1);
    aliasbygroup=true;   % Try to avoid aliasing by group first
    mzsepindex=1;	% Start with largest separation if possible
    nbackup=0;
    while ~any(vecs(:,t))
      % Choose a vector
      % Non-full vector that is not aliased (by mass)
      vbest=find(~valias & sum(vecs,2)<targetspervec);
      if isempty(vbest)
        % No possible choices that aren't aliased
        if nbackup==0
          % First try backing up on a vector that is full already
          fullvecs=find(~valias & sum(vecs,2)==targetspervec);
          if ~isempty(fullvecs)
            % Back up 
            bv=fullvecs(randi(length(fullvecs),1,1));
            [al,mdiff]=aliased(mass(t),mass(vecs(bv,:)));
            if al
              % Wouldn't help
              valias(bv)=true;
              continue;  % Try again
            end
            bts=find(vecs(bv,:));
            bt=bts(randi(length(bts),1,1));
            vecs(bv,bt)=false;
            nbackup=nbackup+1;
            continue;
          end
        end
        if mzsepindex<length(mzminsep)
          mzsepindex=mzsepindex+1;
          valias(:)=false;
          nbackup=0;
          continue;
        end
        if aliasbygroup
          % Switch to aliasing by single vector
          aliasbygroup=false;
          mzsepindex=1;
          valias(:)=false;
          nbackup=0;
          continue;
        end
        % Stuck!
        assert(false);
      end
      if ~isempty(ndiffs)
        % Try to design vectors that increases the number of diffs with other target in the cumulative vector set
        avoid=sum(vecs(vbest,ndiffs(t,:)==min(ndiffs(t,:))),2);  % Number of targets in each vector that we're trying to avoid
        % Only attempt to place in vectors that have the minimum number of targets that we're trying to avoid
        vbest=vbest(avoid==min(avoid));
      end
      v=vbest(randi(length(vbest),1,1));
      % Check for an alias
      % Check across sets of 4 vectors (since we'll mix them by 4's later)
      if aliasbygroup
        vset=floor((v-1)/4)*4+1;
        [al,mdiff]=aliased(mass(t),mass(any(vecs(vset:vset+3,:),1)),mzminsep(mzsepindex));
      else
        [al,mdiff]=aliased(mass(t),mass(vecs(v,:)),mzminsep(mzsepindex));
      end        
      if al & ~valias(v)   % Aliased and not already flagged
        valias(v)=true;
      else
        if debug && (mzsepindex~=1 || ~aliasbygroup || nbackup>0)
          fprintf('[%d:',sum(vecs(:)));
          if aliasbygroup
            fprintf('M');
          else
            fprintf('S');
          end
          fprintf('%.3f',mzminsep(mzsepindex));
          if nbackup>0
            fprintf(',%d',nbackup);
          end
          fprintf(']');
          ndebug=ndebug+1;
          if ndebug>20
            fprintf('\n');
            ndebug=0;
          end
        end
        cnt(t,:,:)=0;
        cnt(t,aliasbygroup+1,mzsepindex)=1;
        vecs(v,t)=true;
        assert(sum(vecs(v,:))<=targetspervec);
        if al
          fprintf('[%d,%.4f]',length(mdiff),min(mdiff));
        end
      end
    end
  end
  m=[];
  for i=1:nvectors
    m(i)=min(diff(sort(mass(vecs(i,:)))));
  end
  fprintf('min delta mass = %.4f\n',min(m));
  ccodes='SM';
  total=squeeze(sum(cnt,1));
  for i=1:size(total,1)
    for j=1:size(total,2)
      fprintf('%c %.3f: %d\n', ccodes(i),mzminsep(j), total(i,j));
    end
  end
end


  