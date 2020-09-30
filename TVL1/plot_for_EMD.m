%%
imageplot(mu1);
hold on;
imageplot(mu2);
hold on;

imageplot(rescale(crop(mu1-mu2,256)),[],'fit');
hold on;
mx=reshape(m([1:end/2]),256,256);
my=reshape(m([end/2+1:end]),256,256);
quiver(my,mx);