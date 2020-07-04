s7sdf=SDF();
s7sdf.read('../../data/HTBCFiles/S7Set.sdf');
s7sdf.set384();
s7sdf.getformulae();
save('../../data/matfiles/S7SDF.mat','s7sdf');


mass=[s7sdf.sdf.MonoisotopicMass];
MASSFUZZ=0.003;
nsame=[];
groupsize=320;
for k=1:10
  sel=randi(length(mass),groupsize,1);
  smass=sort(mass(sel));
distinct=true(length(smass));
for i=1:length(smass)
  for j=i+1:length(smass)
    if smass(j)-smass(i) < 2*MASSFUZZ
      distinct(i,j)=false;
      distinct(j,i)=false;
    else
      break;
    end
  end
end
nsame(k)=mean(sum(~distinct));
end

fprintf('For groups of %d, each molecule has similar mass as %.2f others on average\n', groupsize, mean(nsame));



