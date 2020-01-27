function [map,sigma]=piecewise(x,y,maxrmse,maxsegments)
% Find piecewise linear map from x to y with RMSE <= maxrmse
% Try 1 piece first
minpoints=10;
[x,ia]=unique(x,'rows','sorted');
y=y(ia);
sel=isfinite(x)&isfinite(y);
x=x(sel); y=y(sel);
[fit,stats]=robustfit(x,y);
fprintf('Piecewise over %d points [%f-%f] -> %f\n', length(x), min(x), max(x),stats.robust_s);
assert(fit(2)>0); % Make sure it is monotonic increasing
if stats.robust_s<maxrmse || length(x)<minpoints*2 || maxsegments<2
  map=[min(x),fit(1)+fit(2)*min(x); max(x), fit(1)+fit(2)*max(x)];
  sigma=stats.robust_s;
  return;
end
[x,ord]=sort(x);
y=y(ord);
stotal=[];
for i=minpoints:length(x)-minpoints
  [~,s1]=robustfit(x(1:i),y(1:i));
  [~,s2]=robustfit(x(i+1:end),y(i+1:end));
  stotal(end+1)=sqrt((s1.robust_s^2*i + (length(x)-i)*s2.robust_s^2)/length(x));
  %fprintf('i=%d, s1=%f, s2=%f, stotal=%f\n',i,s1.robust_s, s2.robust_s, stotal(end));
end
[~,i]=min(stotal);
i=i+minpoints-1;
fprintf('Using i=%d\n', i);
[map1,s1]=piecewise(x(1:i),y(1:i),maxrmse,ceil(maxsegments/2));
[map2,s2]=piecewise(x(i+1:end),y(i+1:end),maxrmse,floor(maxsegments/2));
sigma=sqrt((s1^2*i + (length(x)-i)*s2^2)/length(x));
map=[map1(1:end-1,:);(map1(end,:)+map2(1,:))/2;map2(2:end,:)];
fprintf('Split x(1:%d) at %d -> %f, %f -> %f\n', length(x), i, s1, s2, sigma);
end

