% Use hough transform to do a regression
function [m,b]=houghregress(x,y,debug)
if nargin<3
  debug=false;
end
N=500;
im=false(N,N);
xmap=round((x-min(x))/(max(x)-min(x))*(N-1)+1);
ymap=round((y-min(y))/(max(y)-min(y))*(N-1)+1);
im(sub2ind(size(im),ymap,xmap))=true;
[H,T,R]=hough(im,'Theta',[-90:-45,45:89]);
P=houghpeaks(H,2);
theta=T(P(1,2));
rho=R(P(1,1));
mmap=-cosd(theta)/sind(theta);
bmap=rho/sind(theta);
linex=[1,N];
liney=mmap*linex+bmap;


ux=(linex-1)*(max(x)-min(x))/(N-1)+min(x);
uy=(liney-1)*(max(y)-min(y))/(N-1)+min(y);
pfit=polyfit(ux,uy,1);
m=pfit(1);
b=pfit(2);

if debug
  setfig('hough0');clf;
  plot(x,y,'o');
  hold on;
  ax=axis;
  plot(ax(1:2),m*ax(1:2)+b,'r');
  plot(ux,uy,'ro');

  setfig('hough1');clf;
  imshow(im);
  axis xy
  ax=axis;
  hold on;
  plot(linex,liney,'r');

  setfig('hough2');clf;
  imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');25
  xlabel('\theta'), ylabel('\rho');
  axis on, axis normal, hold on;
  plot(T(P(:,2)),R(P(:,1)),'s','color','red');
  keyboard
end

end
