% Well for sample i
% Order in offset 96-well pattern  (A1, A2, B1, B2; C1, C2, D1, D2, ...)
% Then it'll be easier to create V256 96-well plates from these
function s=wellname384(i)
  i96=floor((i-1)/4); 
  w=mod(i-1,4);
  c=floor(i96/8)*2;
  r=mod(i96,8)*2;
  c=c+(w==1 | w==3);
  r=r+(w==2 | w==3);
  s=sprintf('%c%02d',char('A'+r),c+1);
  %  s=sprintf('%c%02d',char('A'+mod(i-1,16)),floor((i-1)/16)+1);
end

