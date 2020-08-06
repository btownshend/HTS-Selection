function [mz,name]=adducts(M,thresh)
if nargin<2
  thresh=.01;
end
% Name, m/z, NIST12 counts
list={'M+3H',M/3 + 1.007276,105,
'M+2H+Na',M/3 + 8.334590,0,
'M+H+2Na',M/3 + 15.7661904,0,
'M+3Na',M/3 + 22.989218,0,
'M+2H',M/2 + 1.007276,18025,
'M+H+NH4',M/2 + 9.520550,0,
'M+H+Na',M/2 + 11.998247,288,
'M+H+K',M/2 + 19.985217,276,
'M+ACN+2H',M/2 + 21.520550,0,
'M+2Na',M/2 + 22.989218,217,
'M+2ACN+2H',M/2 + 42.033823,0,
'M+3ACN+2H',M/2 + 62.547097,0,
'M+H',M + 1.007276,98521,
'M+NH4',M + 18.033823,1161,
'M+Na',M + 22.989218,8679,
'M+CH3OH+H',M + 33.033489,1,
'M+K',M + 38.963158,215,
'M+ACN+H',M + 42.033823,0,
'M+2Na-H',M + 44.971160,81,
'M+IsoProp+H',M + 61.06534,0,
'M+ACN+Na',M + 64.015765,0,
'M+2K-H',M + 76.919040,0,
'M+DMSO+H',M + 79.02122,0,
'M+2ACN+H',M + 83.060370,0,
'M+IsoProp+Na+H',M + 84.05511,0,
'2M+H',2*M + 1.007276,0,
'2M+NH4',2*M + 18.033823,0,
'2M+Na',2*M + 22.989218,0,
'2M+K',2*M + 38.963158,0,
'2M+ACN+H',2*M + 42.033823,0,
'2M+ACN+Na',2*M + 64.015765,0};
cnt=cell2mat(list(:,3));
frac=cnt/sum(cnt);
keep=find(frac>=thresh);
[~,ord]=sort(frac(keep),'descend');
keep=keep(ord);
mz=cell2mat(list(keep,2));
name=list(keep,1);
end

