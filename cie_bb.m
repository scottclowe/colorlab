% Black-body

close all;

% For a pre-existing implementation, see:
% http://www.vendian.org/mncharity/dir3/blackbody/UnstableURLs/tool_pl.txt
% Also of some help is
% http://www.fourmilab.ch/documents/specrend/specrend.c

% Physical constants
% Taken from Wikipedia on 2013-06-26
h = 6.62606957e-34; % Planck's constant, J.s
c = 299792458;      % Speed of light, m/s
k = 1.3806488e-23;  % Boltzmann's constant, J/K

% sRGB primaries
%            x      y      Y
primaries = [0.6400 0.3300 0.2126;   % Red
             0.3000 0.6000 0.7153;   % Green
             0.1500 0.0600 0.0721;]; % Blue
wp        = [0.3127 0.3290 1.0000 ]; % White

% Color matching functions taken from http://www.cvrl.org/cmfs.htm
cie_cmf_fname = 'lin2012xyz2e_fine_7sf.csv';
cie_cmf = csvread(cie_cmf_fname);

figure; hold on;
plot(cie_cmf(:,1),cie_cmf(:,2),'r');
plot(cie_cmf(:,1),cie_cmf(:,3),'g');
plot(cie_cmf(:,1),cie_cmf(:,4),'b');
xlabel('Wavelength (nm)');
title('Color matching functions');

% Extract wavelengths from CMF file
lambda = cie_cmf(:,1); % Column vector
lambda = lambda*1e-9;  % nm -> m

% Have to weight for x-spacing of bars
dlambda = (lambda(2)-lambda(1)) * 1e9;

% Do the whitepoint D65 temperature
T65 = 6503.6;
I65 = 2*pi*h*c^2 ./ bsxfun(@times, lambda.^5, (exp( h*c./k *(1./(lambda*T65)) ) -1 ) );
XYZ65 = I65' * cie_cmf(:,2:4) * dlambda;
xyz65 = XYZ65/sum(XYZ65);

% T = 0:1000:10000; % Temperatures in K, Row vector
% Make our set of T distributed evenly in reciprocal space (aka mired)
nT = 100;
T1 = 1000;
T2 = 30000;
mired = linspace(1e6/T1,1e6/T2,nT);
T = 1e6./mired;
% Also include some low temperatures
T = [100:10:990 T];

% Generate spectral intensities using Planck's law
I = 2*pi*h*c^2 ./ bsxfun(@times, lambda.^5, (exp( h*c./k *(1./(lambda*T)) ) -1 ) );

% Normalise I so max is 1
% I = bsxfun(@rdivide,I,max(I,[],1));

% Normalise I so that it integrates to be 1
% I = bsxfun(@rdivide,I,sum(I,1)*dlambda);
% % I = I*100;

% Integrate with rectangle approximation
XYZ = I' * cie_cmf(:,2:4) * dlambda;

% Normalise Y so the brightest point is at Y=1?
% XYZ = bsxfun(@rdivide,XYZ,max(XYZ(:,2));

% Normalise XYZ so D65 has Y=1
% Do the whitepoint D65 temperature
T65 = 6503.6;
I65 = 2*pi*h*c^2 ./ bsxfun(@times, lambda.^5, (exp( h*c./k *(1./(lambda*T65)) ) -1 ) );
XYZ65 = I65' * cie_cmf(:,2:4) * dlambda;
xyz65 = XYZ65/sum(XYZ65);
% XYZ = XYZ/XYZ65(2);
XYZ = bsxfun(@rdivide,XYZ,max(XYZ(:,2),XYZ65(2)));

figure;
plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'-o');
xlabel('X');
ylabel('Y');
zlabel('Z');

xyz = bsxfun(@rdivide,XYZ,sum(XYZ,2));

figure; hold on;
plot(xyz(:,1),xyz(:,2),'-o');
plot(primaries([1:3 1],1),primaries([1:3 1],2),'k-');
plot(primaries(1,1),primaries(1,2),'ro');
plot(primaries(2,1),primaries(2,2),'go');
plot(primaries(3,1),primaries(3,2),'bo');
plot(wp(1),wp(2),'kx');
xlabel('x');
ylabel('y');
xlim([0 1]);
ylim([0 1]);
set(gca,'Color',[0.4663 0.4663 0.4663]);

figure;
plot3(xyz(:,1),xyz(:,2),XYZ(:,2),'-o');
xlabel('x');
ylabel('y');
zlabel('Y');

% Normalise intensity by finding brightest, most saturated point with x,y
Y = find_gamut_xy_maxY(xyz(:,1),xyz(:,2),primaries');

xyY = xyz;
xyY(:,3) = Y;

fprintf('%d out of %d points out of gamut\n',sum(isnan(Y)),length(Y));

xyY = xyY(~isnan(Y),:);

xyz2 = [xyY(:,1) xyY(:,2) 1-xyY(:,1)-xyY(:,2)];
XYZ2 = bsxfun(@times,xyz2,xyY(:,3)./xyY(:,2));
% XYZ2 = applycform([x y repmat(0.2,size(y))],makecform('xyl2xyz'));

% figure;
% plot3(xyY(:,1),xyY(:,2),xyY(:,3),'-o');
% xlabel('x');
% ylabel('y');
% zlabel('Y');
% 
% figure;
% plot3(XYZ2(:,1),XYZ2(:,2),XYZ2(:,3),'-o');
% xlabel('X');
% ylabel('Y');
% zlabel('Z');


rgb = applycform(XYZ2,makecform('xyz2srgb'));

clrs = repmat(rgb,[1 1 20]);
clrs = permute(clrs,[1 3 2]);
figure;
imagesc(clrs);
axis xy;

Lab = applycform(XYZ2,makecform('xyz2lab'));
% cform = makecform('xyl2lab');
% Lab = applycform([x y Y],cform);

figure;
plot3(Lab(:,2),Lab(:,3),Lab(:,1),'-o')
xlabel('a*')
ylabel('b*')
zlabel('L*')
set(gca,'Color',[0.4663 0.4663 0.4663]);

% g = fetch_cielchab_gamut;
% % Get a mesh version of the gamut
% if ~isfield(g,'lchmesh')
%     g.lchmesh = make_gamut_mesh(g);
% end
% Lm = g.lchmesh.Lgrid;
% cm = g.lchmesh.cgrid;
% hm = g.lchmesh.hgrid/180*pi;
% am = cm.*cos(hm);
% bm = cm.*sin(hm);
% 
% cform = makecform('lab2srgb');
% rgbm = applycform([Lm(:) am(:) bm(:)], cform);
% chart = 1:size(rgbm,1);
% chart = reshape(chart,size(Lm));
% 
% figure;
% mesh(am,bm,Lm,chart,'facealpha',0.5,'edgealpha',0.5);
% colormap(rgbm);
% hold on;
% plot3(Lab(:,2), Lab(:,3), Lab(:,1), 'k-');
% set(gca,'Color',[0.4663 0.4663 0.4663]);
% xlabel('a*');
% ylabel('b*');
% zlabel('L*');
% view(0,90);