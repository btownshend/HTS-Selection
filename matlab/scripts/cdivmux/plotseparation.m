% Plot separation of multiple hits using vecs
function plotseparation(vecs,usemax)
  if nargin<2
    usemax=false;
  end
   
  nboot=100;
  ti=sprintf('Separation V%d',mean(sum(vecs,2)));
  setfig(ti);clf;
  leg={};h=[];
  for i=1:10
    for b=1:nboot
      p=randperm(size(vecs,2));
      hits=p(1:i);
      pos=true(size(vecs));  % Possible hits
      for k=1:size(vecs,1)
        if k>1
          pos(k,:)=pos(k-1,:);
        end
        if ~any(vecs(k,hits))
          pos(k,vecs(k,:))=false;
        end
      end
      npos(:,b)=sum(pos,2);
    end
    if usemax
      v=max(npos,[],2);
    else
      v=mean(npos,2);
    end
    h(i)=plot(1:size(vecs,1),v);
    leg{end+1}=sprintf('%d hits',i);
    hold on;
    enough=find(v<=i*2,1);
    if ~isempty(enough)
      plot(enough,v(enough),'o','Color',get(h(i),'Color'));
    end
  end
  legend(h,leg);
  xlabel('Num Vectors');
  ylabel('Num possible targets');
  title(ti);
  ax=axis;
  ax(1)=20;
  ax(4)=200;
  axis(ax);
end

