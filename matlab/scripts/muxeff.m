% Compute efficiency of mux
dlen=[8,10,8];
maxhits=20;

nposs=[];
for ndiag=0:3
  for nhits=1:maxhits
    nposs(ndiag+1,nhits)=1-prod(1./dlen(1:ndiag))+prod(1-((dlen(1:ndiag)-1)./dlen(1:ndiag)).^nhits)*(1-(7/8)^nhits)*(1-(11/12)^nhits)*(1-(9/10)^nhits)*960/nhits;
  end
end

ti=sprintf('MUX efficiency - Diag len %d,%d,%d',dlen)
setfig(ti);clf;
h=plot(nposs',':');
xlabel('Hits');
ylabel('Avg. possible targets/Actual');


% And MC method
id=(1:960)';
plate=floor((id-1)/80)+1;
col=floor(mod((id-1),80)/8)+1;
row=mod(id-1,8)+1;
d1=mod(row+plate,dlen(1))+1;
d2=mod(col+plate,dlen(2))+1;
d3=mod(row+col,dlen(3))+1;
alltgts=[plate,col,row,d1,d2,d3];

% Generate some random hits
nmc=900;   % Less than 960 to get unique ones
row=randi(8,nmc,1);
col=randi(10,nmc,1);
plate=randi(12,nmc,1);
d1=mod(row+plate,dlen(1))+1;
d2=mod(col+plate,dlen(2))+1;
d3=mod(row+col,dlen(3))+1;
id=row+(col-1)*8+(plate-1)*80;
meas=[plate,col,row,d1,d2,d3];

% Make sure there are no repeats within maxhits of each other
[id,sel]=unique(id,'stable');
meas=meas(sel,:);

cnt=[];
for ndiag=0:3
  for nhits=1:maxhits
    for i=1:size(meas,1)-maxhits+1
      pos=true(960,1);
      for k=1:3+ndiag
        pos(~ismember(alltgts(:,k),meas(i:i+nhits-1,k)))=false;
      end
      cnt(ndiag+1,nhits,i)=sum(pos);
    end
  end
end
nposs=mean(cnt,3);
for nhits=1:maxhits
  nposs(:,nhits)=nposs(:,nhits)/nhits;
end

hold on;
for i=1:length(h)
  plot(nposs(i,:)','o-','Color',get(h(i),'Color'));
end

title(ti);
nms=arrayfun(@(z) sprintf('NDiag=%d (est)',z),0:3,'UniformOutput',false);
nms=[nms,arrayfun(@(z) sprintf('NDiag=%d (MC)',z),0:3,'UniformOutput',false)];
legend(nms,'location','best');

set(gca,'YScale','log');
ax=axis;
ax(3)=1;ax(4)=20;
axis(ax);
logticks(1,1);
