plot3(xA,yA,zA,'-b',xB,yB,zB,'-r');
title('Plot Title')
xlabel('X data')
ylabel('Y data')
zlabel('Z data')
legend('Data A','Data B')
grid on
az = 0;
el = 90;
view([az,el])
degStep = 5;
detlaT = 0.1;
fCount = 71;
f = getframe;
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,fCount) = 0;
k = 1;
% spin 45°
for i = 0:-degStep:-45
  az = i;
  ([az,el])
  f = getframe;
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
  k = k + 1;
end
% tilt down
for i = 90:-degStep:15
  el = i;
  view([az,el])
  f = getframe;
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
  k = k + 1;
end
% spin left
for i = az:-degStep:-90
  az = i;
  view([az,el])
  f = getframe;
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
  k = k + 1;
end
% spin right
for i = az:degStep:0
  az = i;
  view([az,el])
  f = getframe;
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
  k = k + 1;
end
% tilt up to original
for i = el:degStep:90
  el = i;
  view([az,el])
  f = getframe;
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
  k = k + 1;
end
imwrite(im,map,'Animation.gif','DelayTime',detlaT,'LoopCount',inf)