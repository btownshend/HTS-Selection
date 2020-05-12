% Build a logic matrix size (nvectors,ntargets) such that each row has targetspervec true values
% Total number of trues in each column should be equal within 1
function [vecs,vecspertarget]=buildvecs(nvectors,ntargets,targetspervec)
  vecs=false(nvectors,ntargets);
  if length(targetspervec)==1
    targetspervec=repmat(targetspervec,nvectors,1);
  end
  assert(length(targetspervec)==nvectors);
    
  % Fill in row-by-row
  for v=1:nvectors
    while sum(vecs(v,:))<targetspervec(v)
      % Choose one of the columns that are most needed at this point
      % Has the side effect that vectors 1:k will give uniform coverage of the targets (+-1)
      tbest=find(sum(vecs,1)==min(sum(vecs,1)));
      t=tbest(randi(length(tbest),1,1));
      if vecs(v,t)
        % fprintf('.');   % Already taken
        ;
      else
        vecs(v,t)=true;
        %fprintf('+');
      end
    end
    %fprintf('\n');
  end
  
  vecspertarget=sum(targetspervec)/ntargets;
  assert(all(sum(vecs)<=ceil(vecspertarget) & sum(vecs)>=floor(vecspertarget)));
  assert(all(sum(vecs,2)==targetspervec));
end
