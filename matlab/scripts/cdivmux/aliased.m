% Check if a given mass may be aliased with any of the others
function [al,diff]=aliased(mass,others,mzminsep)
% Possible adducts (with z=1)
  ADDUCTS=[1.007276,18.033823,22.989218,33.033489,38.963158,42.033823,44.971160,61.065340,64.015765,76.919040,79.021220,83.060370,84.055110];
  if nargin<3
    mzminsep=0.01;
  end
  if isempty(others)
    al=false;
    diff=inf;
  else
    tgts=[];
    diff=ones(size(others))*1e10;
    for i=1:length(ADDUCTS)
      diff=min(diff,abs((mass+ADDUCTS(1))-(others+ADDUCTS(i))));
      diff=min(diff,abs((mass+ADDUCTS(i))-(others+ADDUCTS(1))));
    end
    al=any(diff<mzminsep);
  end
end
