function cmap = cie_bluewhitered_cmap2(n, func, debug)

if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1);
end
if nargin<2
    func = [];
end
if nargin<3
    debug = false;
end

% CIELCH      [  L*    c    h]
  lchblue   = [ 44.25,  91.70, 292]; % 289
  lchred    = [ 44.25,  91.70,  41];
  wp        = [ 99.25,   0   ,   0];

% CIELab    [  L*   a*   b*]
labblue   = [lchblue(1), lchblue(2)*cosd(lchblue(3)), lchblue(2)*sind(lchblue(3))];
labred    = [lchred(1) ,  lchred(2)*cosd(lchred(3)) ,  lchred(2)*sind(lchred(3)) ];
labwp     = [wp(1)     ,      wp(2)*cosd(wp(3))     ,      wp(2)*sind(wp(3))     ];

neach = ceil(n/2);

L1 = linspace(labblue(1), labwp(1), neach);
a1 = linspace(labblue(2), labwp(2), neach);
b1 = linspace(labblue(3), labwp(3), neach);

L2 = linspace(labwp(1), labred(1), neach);
a2 = linspace(labwp(2), labred(2), neach);
b2 = linspace(labwp(3), labred(3), neach);

Lab1 = [L1' a1' b1'];
Lab2 = [L2' a2' b2'];

switch mod(n,2)
    case 0
        % n is even
        Lab1 = Lab1(1:end-1,:);
        Lab2 = Lab2(2:end  ,:);
    case 1
        % n is odd
        Lab1 = Lab1(1:end-1,:);
%         Lab2 = Lab2;
end

Lab  = [Lab1;Lab2];

if ~isempty(func)
    cmap = func(Lab);
elseif license('checkout','image_toolbox')
    % If using ImageProcessingToolbox
    cform = makecform('lab2srgb');
    cmap = applycform(Lab, cform);
elseif exist('colorspace','file')
    % Use colorspace
%     warning('LABWHEEL:NoIPToolbox:UseColorspace',...
%         ['Could not checkout the Image Processing Toolbox. ' ...
%          'Using colorspace function.']);
    cmap = colorspace('Lab->RGB',Lab);
else
    % Use colorspace
    warning('LABWHEEL:NoIPToolbox:NoColorspace',...
        ['Could not checkout the Image Processing Toolbox. ' ...
         'Colorspace function not present either.\n' ...
         'You need one of the two to run this function.']);
     if exist('suggestFEXpackage','file')
        suggestFEXpackage(28790,'You may wish to download the colorspace package.');
     end
end


if debug
    
g = fetch_cielchab_gamut('srgb');

% Finding the right values
% h = lchred(3); %292;
% li = g.lch(:,3)==h;
% gh = g.lch(li,:);
% 
% figure; hold on; set(gca,'Color',[.48 .48 .48]);
% scatter(gh(:,2),gh(:,1),20,g.rgb(li,:));
% M = (100-gh(end-6,1))/(0-gh(end-6,2));
% C = 100; Y = 100:-1:40; X = (Y-C)/M;
% plot(X,Y,'k-');


% Checking right values

% li = g.lch(:,3)==lchblue(3);
% gh = g.lch(li,:);
% figure; 
% hold on;
% subplot(1,2,1);
% % scatter(gh(:,2),gh(:,1),20,g.rgb(li,:));
% plot(gh(:,2),gh(:,1),'b-');
% set(gca,'Color',[.48 .48 .48]);
% plot([0 lchblue(2)],[100 lchblue(1)],'o-k');
% box on;
% 
% li = g.lch(:,3)==lchred(3);
% gh = g.lch(li,:);
% subplot(1,2,2);
% % scatter(gh(:,2),gh(:,1),20,g.rgb(li,:));
% plot(gh(:,2),gh(:,1),'r-');
% set(gca,'Color',[.48 .48 .48]);
% hold on;
% plot([0 lchred(2)],[100 lchred(1)],'o-k');
% box on;


ghb = g.lch(g.lch(:,3)==lchblue(3),:);
ghr = g.lch(g.lch(:,3)==lchred(3),:);

figure; set(gca,'Color',[.48 .48 .48]); hold on; box on;
plot(ghb(:,2),ghb(:,1),'b-');
plot([labwp(2) lchblue(2)],[labwp(1) lchblue(1)],'k-');
plot(ghr(:,2),ghr(:,1),'r-');
plot([labwp(2) lchred(2)],[labwp(1) lchred(1)],'ko');


% Plot the colormap
img = repmat(cmap,[1 1 20]);
img = permute(img,[1 3 2]);
figure;
imagesc(img);
axis xy;
% figure;
% imagesc(img(1:4:end,:,:));
% axis xy;
end