% Groups
ngroups=8;
vectorsperfinal=3;   % How many per-group vectors go into each final vector; if we want ~20uM, then this should be ~250/80
nfinalvectors=96*4;   % Make it such that every group vector is used equally

% Per-group stats
nplatespergroup=8;
ntargetsperplate=80;
ntargets=nplatespergroup*ntargetsperplate;
nvectors=96;
targetspervec=80;
nvecspertarget=nvectors*targetspervec/ntargets;

if ~exist('vecs','var')
  if ngroups==1
    vecs=buildvecs(nvectors,ntargets,targetspervec);
  else
    gvecs=[];
    fprintf('Building %d group vectors...',nvectors);
    for i=1:ngroups
      gvecs(i,:,:)=buildvecs(nvectors,ntargets,targetspervec);
    end
    fprintf('done\n');
    
    % Build final vecs using picking from gvecs
    fprintf('Building %d final picks...',nfinalvectors);
    finalpick=buildvecs(nfinalvectors,ngroups,vectorsperfinal);

    % Combine into a single overall vecs matrix
    vecs=false(nfinalvectors,ntargets*ngroups);
    nextvector=ones(ngroups,1);  % Next vector to use from each group
    for i=1:nfinalvectors
      for j=1:ngroups
        if finalpick(i,j)
          sel=(1:ntargets)+(j-1)*ntargets;
          vecs(i,sel)=squeeze(any(gvecs(j,nextvector(j),:),2))';
          nextvector(j)=mod(nextvector(j),nvectors)+1;
          assert(sum(vecs(i,sel),2)==targetspervec);
        end
      end
    end
    fprintf('done\n');
    fprintf('Each vector has %d/%d targets\n', mean(sum(vecs,2)),size(vecs,2));
    fprintf('Each target is in %.2f/%d vectors\n', mean(sum(vecs,1)),size(vecs,1));
    assert(max(sum(vecs,1))-min(sum(vecs,1))<=1);
    assert(all(abs(sum(vecs,2)-targetspervec*vectorsperfinal)<=1));
  end
end

% Run some tests
mainti=sprintf('randvecs: ngroups=%d, targetspervec=%d, vectorsperfinal=%d', ngroups, targetspervec,vectorsperfinal);
setfig(mainti);clf;
h=tiledlayout('flow');
title(h,mainti);

nmc=100;   % Number of monte carlo simulations
step=47;
while size(vecs,1)/step > 8
  step=step+47;
end
for nvecs=step:step:size(vecs,1)
  npos=[];
  ncert=[];
  for nhits=1:3*(ngroups+4)   % Number of targets/aptamers
    phot=1-(1-targetspervec/ntargets)^nhits;
    for rep=1:nmc
      [pos,cert]=simvecs(vecs,nhits,nvecs);
      npos(nhits,rep)=sum(pos);
      ncert(nhits,rep)=sum(cert);
    end
    
    fprintf('nhits=%d, phot=%.2f, npos=%.1f+=%.1f, ncert=%.1f+-%.1f\n', ...
            nhits, phot, mean(npos(nhits,:)),std(npos(nhits,:)),mean(ncert(nhits,:)),std(ncert(nhits,:)));
  end

  ti=sprintf('nvectors=%d', nvecs);
  nexttile;
  boxplot(npos');
  set(gca,'YScale','log');
  hold on;
  plot(1:size(npos,1),1:size(npos,1),':g');
  ax=axis(); ax(3:4)=[1,100];
  axis(ax);
  ylabel('Num of possible targets');

  yyaxis right
  plot(mean(ncert,2));
  set(gca,'YScale','log');
  ylabel('Num of certain targets');
  axis(ax);
  
  xlabel('Number of simultaneous hits');
  set(gca,'XTickLabelRotation',45);
  title(ti)
  pause(0.1);
end

