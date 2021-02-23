% Check if plate 565 was inverted

expect=find(strcmp(v256.names,'565A02')):find(strcmp(v256.names,'565H11'));  % Compounds in plate 565
actual=find(strcmp(v256.names,'568A02')):find(strcmp(v256.names,'568H11'));  % Compounds in plate 565
                                                                             %actual=actual(end:-1:1);
adduct=1;
tests=[];
for i=1:length(cc)
  c=expect(i);
  c2=actual(i);
  %c2=c;
  fprintf('Looking for %s in vectors that contain %s:\n', v256.names{c},v256.names{c2});
  for j=1:length(v256.samples)
    if v256.contains(c2,j);   % Check the vectors that are supposed to contain the flipped compound
      mh=v256.multihits{c,adduct,j};
      fprintf(' %s(%d): %s\n', v256.samples{j},j, strjoin(arrayfun(@(z) v256.reffeatures(j).features(z).tostring(),mh,'Unif',false),'; '));
      tests(end+1)=~isempty(mh);
    end
  end
end
fprintf('Passed %d/%d tests\n', sum(tests), length(tests));
