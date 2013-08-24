
use_uplab = false;

g = fetch_cielchab_gamut('srgb', [], [], use_uplab);

% Start values
npoints = 100;
expnt = 2;
L_off = 0;
L0 = 5;
Lmax = 95;
% c0 = 0;
typ = 'pow';

L = linspace(L0,Lmax,npoints);
Lmid = (L0+Lmax)/2;

% hue1_range = 0:359;
% hue2_range = -360:720;

hue1_range =    0:2:359;
hue2_range = -360:2:720;

all_maxc = nan(length(hue1_range),length(hue2_range));

for ih1=1:length(hue1_range)
    for ih2=1:length(hue2_range);
        h1 = hue1_range(ih1);
        h2 = hue2_range(ih2);
        h_per_L = (h2-h1)/(Lmax-L0);

        h = h1 + h_per_L * ((L-L0) + L_off - (L-Lmid).^2 * L_off/(Lmax-Lmid)^2);
        h = mod(h,360);
        
        switch typ
            case 'sin'
                c = sin(pi* (L-L0)/(Lmax-L0) ).^expnt;
            case 'pow'
                c = 1 - abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(Lmax-Lmid,Lmid-L0))).^expnt) / abs(Lmid.^expnt);
            otherwise
                error('Unfamiliar type');
        end

        Lch = [L' c' h'];

        % Check for points out of gamut
        [TF,P2] = isingamut(Lch,g,'Lch');

        P2C = P2(:,2);
        maxc = min(P2C./c');
        
        all_maxc(ih1,ih2) = maxc;
    end
end

%%

figure;
imagesc(hue2_range, hue1_range, all_maxc);
colormap(cie_hot);
colorbar;
ylabel(sprintf('Start hue at L=%.1f',L0));
xlabel(sprintf('End hue at L=%.1f',Lmax));

%%
return;

%%

% h1 = 275;
% h2 = 360;

h1 = 290;
h2 = 380;

h_per_L = (h2-h1)/(Lmax-L0);

h = h1 + h_per_L * ((L-L0) + L_off - (L-Lmid).^2 * L_off/(Lmax-Lmid)^2);
h = mod(h,360);

switch typ
    case 'sin'
        c = sin(pi* (L-L0)/(Lmax-L0) ).^expnt;
    case 'pow'
        c = 1 - abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(Lmax-Lmid,Lmid-L0))).^expnt) / abs(Lmid.^expnt);
    otherwise
        error('Unfamiliar type');
end

Lch = [L' c' h'];

% Check for points out of gamut
[TF,P2] = isingamut(Lch,g,'Lch');

P2C = P2(:,2);
maxc = min(P2C./c');

c = c * maxc;


Lch = [L' c' h'];
a = c.*cosd(h);
b = c.*sind(h);
Lab = [L' a' b'];

rgb = hard_lab2rgb(Lab, use_uplab);

img = repmat(rgb,[1 1 20]);
img = permute(img,[1 3 2]);
figure;
imagesc(img);

plot_labcurve_rgbgamut(Lab, use_uplab);

%%
% Redo with maxc3
c = maxc * c;

% Go from cylindrical Lch co-ords to cartesian Lab
a = c.*cos(h/360*(2*pi));
b = c.*sin(h/360*(2*pi));

% Turn Lab into sRGB values
Lab = [L' a' b'];
rgb = soft_lab2rgb(Lab, use_uplab);

% Check for points out of gamut
[TF,P2] = isingamut(Lab,g,'Lab');
fprintf('Method #%d now has %d of %d in gamut.\n',method,sum(TF),length(TF));

P2C = sqrt(P2(:,2).^2+P2(:,3).^2);


% First and last colors may be out due to floating point error
% others should not be
% Lab(~TF,:) = P2(~TF,:);

% Plot map and its derivation
img = repmat(rgb,[1 1 20]);
img = permute(img,[1 3 2]);

figure;
imagesc(img);
axis xy;
title(sprintf('expnt = %.3f; h1 = %.3f; h/L = %.3f; L.off = %.3f',...
    expnt, h1, h_per_L, L_off));

figure;
hold on;
plot(c1, L, 'Color', [0 .8 0]);
plot(c , L, 'k');
plot(sqrt(P2(:,2).^2+P2(:,3).^2), P2(:,1), 'r');
title(sprintf('expnt = %.3f; h1 = %.3f; h/L = %.3f; L.off = %.3f',...
    expnt, h1, h_per_L, L_off));
xlabel('chroma');
ylabel('Lightness');

figure;
hold on;
plot(L,c,'k');
plot(L,h,'r');
title(sprintf('expnt = %.3f; h1 = %.3f; h/L = %.3f; L.off = %.3f',...
    expnt, h1, h_per_L, L_off));
xlabel('Lightness');
legend('chroma','hue','Location','EastOutside')

figure;
scatter3(Lab(:,2), Lab(:,3), Lab(:,1), 20, rgb, 'filled');
set(gca,'Color',[0.4663 0.4663 0.4663]);
xlabel('a');
ylabel('b');
zlabel('L');

figure;
hold on;
plot3(Lab(:,2), Lab(:,3), Lab(:,1), 'bx')
plot3(P2(1:end-1,2) , P2(1:end-1,3) , P2(1:end-1,1) , 'r-') % Last point is grey so P has hue 0
xlabel('a');
ylabel('b');
zlabel('L');
view(0,90);


res = 17;
figure('Position',[150 150 550 550]);
hold on;
plot3(Lab(:,2), Lab(:,3), Lab(:,1), 'k-');
% Lab_s = g.lab(g.lab(:,3)>=0,:);
% rgb_s = g.rgb(g.lab(:,3)>=0,:);
Lab_s = g.lab;
rgb_s = g.rgb;
scatter3(...
    Lab_s([1:res:end-1 end],2), ...
    Lab_s([1:res:end-1 end],3), ...
    Lab_s([1:res:end-1 end],1), ...
    20, ...
    rgb_s([1:res:end-1 end],:), ...
    'filled');
set(gca,'Color',[0.4663 0.4663 0.4663]);
set(gca,'XLim',[-150 150],'YLim',[-150 150],'ZLim',[0 100]);
xlabel('a*')
ylabel('b*')
zlabel('L*')

% % Get a mesh version of the gamut
% if ~isfield(g,'lchmesh')
%     g.lchmesh = make_gamut_mesh(g);
% end
% Lm = g.lchmesh.Lgrid(g.lchmesh.hgrid(:,1)<=180,:);
% cm = g.lchmesh.cgrid(g.lchmesh.hgrid(:,1)<=180,:);
% hm = g.lchmesh.hgrid(g.lchmesh.hgrid(:,1)<=180,:)/180*pi;
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
% xlabel('a');
% ylabel('b');
% zlabel('L');

plot_labcurve_rgbgamut(Lab);