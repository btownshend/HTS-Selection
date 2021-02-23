nout=[];
for i=1:length(res)
  if res(i).adduct==3
    for j=1:length(res(i).samples)
      s1=res(i).samples(j);
      s2=res2(i).samples(j);
      kisos=find(cellfun(@(z) strcmp(z(end-4:end),'[22]K'),{s2.isotopes.name}));
      for k=1:length(kisos)
        r2=s2.isotopes(kisos(k));
        mdiff2=r2.mass-min([s2.isotopes.mass]);
        mdiff1=[s1.isotopes.mass]-min([s1.isotopes.mass]);
        [~,nonkiso]=min(abs(mdiff2-mdiff1));
        assert(length(nonkiso)==1);
        if abs(mdiff2-mdiff1)>0.5
          continue;
        end
        r1=res(i).samples(j).isotopes(nonkiso);
        nout(end+1,:)=[r1.outlier,r2.outlier];
        if nout(end,1)~=nout(end,2)
          fprintf('%d,%d %d=%.3f/%.3f vs %d=%.3f/%.3f\n',i,j,nonkiso,r1.obs,r1.abundance,kisos(k),r2.obs,r2.abundance);
          keyboard;
        end
      end
    end
  end
end