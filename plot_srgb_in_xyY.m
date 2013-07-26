
rgb2xyz = [ ...
 0.4124564  0.3575761  0.1804375;
 0.2126729  0.7151522  0.0721750;
 0.0193339  0.1191920  0.9503041];

N = 51;
Ntot = 6*N.^2;
rgb = nan(Ntot,3);
x = linspace(0,1,N);
[X,Y] = meshgrid(x,x);
rgb0 = [X(:) Y(:) zeros(N.^2,1)];
rgb1 = [X(:) Y(:)  ones(N.^2,1)];
for i=0:2
    j = i*2;
    rgb(j*N.^2+1:(j+1)*N.^2, :) = circshift(rgb0,[0 i]);
    j = i*2+1;
    rgb(j*N.^2+1:(j+1)*N.^2, :) = circshift(rgb1,[0 i]);
end

xyz = rgb * rgb2xyz';

figure;
scatter3(xyz(:,1),xyz(:,2),xyz(:,3),10,rgb,'filled');
set(gca,'Color',[0.4663 0.4663 0.4663]);
xlabel('X');
ylabel('Y');
zlabel('Z');

xyY = bsxfun(@rdivide,xyz,sum(xyz,2));
xyY(:,3) = xyz(:,2);
figure;
scatter3(xyY(:,1),xyY(:,2),xyY(:,3),10,rgb,'filled');
set(gca,'Color',[0.4663 0.4663 0.4663]);
xlabel('x');
ylabel('y');
zlabel('Y');
xlim([0 0.7]);
ylim([0 0.7]);
zlim([0 1]);