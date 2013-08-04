
dbg = 1;
if dbg;
    close all;
    clear all;
    dbg = true;
end;

method = 104;

use_uplab = false;

g = fetch_cielchab_gamut();

% Start values
expnt = 2;
h0 = 0;
h_per_L = 1;
L_off = 0;
maxc = [];
L0 = 0;
Lmax = 100;
c0 = 0;

% Best: method 1
% Would be great if method 5 worked correctly
% Second best: method 8

% Updated:
% Best is one of 25 9 10 11 12

switch method
    case 0
        % Method 1:
        % Caution: goes outside gamut
        %          [  L  c  h  ]
        % Red max chroma
        maxc_red = [54.5 106 40];
        % Yellow max chroma
        maxc_yel = [95.75 94 106];
        % Picked by using cietools to browse with chroma setting

        h_per_L = (maxc_yel(3)-maxc_red(3)) / (maxc_yel(1)-maxc_red(1));
        h0 = maxc_red(3)-h_per_L*maxc_red(1);
        
%         % % Outcome Programmatically:
%         expnt = 2;
%         h0 = -47.2;
%         h_per_L = 1.6;
%         L_off = 0;
%         maxc = 75.85;

    % Method 2: Pick values and play with them
    case 1
        % Set#1
        % Still best at the moment
        expnt = 2;
        h0 = 0;
        h_per_L = 1.15;
        L_off = 0;
        maxc = 70.3;
        
        hmax = h0+100*h_per_L;
        
    case 2
        % Set#2
        expnt = 2;
        h0 = 0;
        h_per_L = 1;
        L_off = 20;
        maxc = 59;
        
        hmax = h0+100*h_per_L;
        
    case 3
        % Set#3
        expnt = 2;
        h0 = 0;
        h_per_L = .9;
        L_off = 25;
        maxc = 60.5;
        
        hmax = h0+100*h_per_L;
    
    case 4
        expnt = 2;
        h0 = -47.2;
        h_per_L = 1.6;
        L_off = 25;
        maxc = 59.3;
        
        hmax = h0+100*h_per_L;

    case 5
        % Would be great if it didn't exceed bounds
        expnt = 2;
%         h0 = 0;
        h0 = 5;
        hmax = 106;
        L_off = 0;
        maxc = []; %75.6;
        
        h_per_L = (hmax-h0)/100;
        
    case 6
        expnt = 2;
        h0 = 0;
        hmax = 106;
        L_off = 0;
        maxc = 68;
        
        h_per_L = (hmax-h0)/100;
        
    case 7
        expnt = 2;
        h0 = 0;
        hmax = 120;
        L_off = -5;
        maxc = 70;
        
        h_per_L = (hmax-h0)/100;
        
    case 8
        expnt = 2.5;
        h0 = 0;
        hmax = 90;
        L_off = 25;
        maxc = 60.5241;
        
        h_per_L = (hmax-h0)/100;
        
    case 21
        % Based on Set#1
        expnt = 2;
        h0 = 0;
        L_off = 0;
        L0   = 2;
        Lmax = 98;
        h_per_L = 1.15;
        maxc = [];
        
        hmax = h0+(Lmax-L0)*h_per_L;

    case 25
        % Based on Set#5
        expnt = 2;
        h0 = 6.5;
        hmax = 106;
        L_off = 0;
        L0   = 2;
        Lmax = 98;
        maxc = 73.9; % 75.6 -> 74.25 -> 73.9
        
        h_per_L = (hmax-h0)/100;
        
    case 28
        % Based on Set#8
        expnt = 2.5;
        h0 = 0;
        hmax = 90;
        L_off = 25;
        L0   = 2;
        Lmax = 98;
        maxc = [];
        
        h_per_L = (hmax-h0)/100;
        
    case 9
        expnt = 2.2;
        h0 = 15;
        hmax = 100;
        L_off = 0;
        L0   = 3;
        Lmax = 97;
        maxc = 73; %[];
        
        h_per_L = (hmax-h0)/100;
        
    case 10
        expnt = 2.5;
        h0 = 15;
        hmax = 100;
        L_off = 5;
        L0   = 3;
        Lmax = 97;
        maxc = 67.8; %[];
        
        h_per_L = (hmax-h0)/100;
        
    case 11
        expnt = 2.3;
        h0 = 10;
        hmax = 100;
        L_off = 5;
        L0   = 3;
        Lmax = 97;
        maxc = 71; %[];
        
        h_per_L = (hmax-h0)/100;
        
    case 12
        expnt = 2.3;
        h0 = 15;
        hmax = 100;
        L_off = 5;
        L0   = 3;
        Lmax = 97;
        maxc = 67.7; %[];
        
        h_per_L = (hmax-h0)/100;
        
    % BLUES (cool)
    case 101
        % Based on Set#25
        expnt = 2;
        h0 = 380;
        hmax = 180;
        L_off = 0;
        L0   = 2;
        Lmax = 98;
        maxc = [];
        
        h_per_L = (hmax-h0)/100;
        
    case 102
        % Based on Set#101
        expnt = 2;
        h0 = 380;
        hmax = 215;
        L_off = 0;
        L0   = 5;
        Lmax = 90;
        maxc = [];
        
        h_per_L = (hmax-h0)/100;
        
    case 103
        % Based on Set#101
        expnt = 2;
        h0 = 375;
        hmax = 172;
        L_off = 0;
        L0   = 5;
        Lmax = 95;
        maxc = 59.7; %[];
        
        h_per_L = (hmax-h0)/100;
        
    case 104
        % Based on Set#103
        expnt = 2;
        h0 = 340;
        hmax = 172;
        L_off = -15;
        L0   = 3;
        Lmax = 85;
        maxc = 70; %[];
        
        h_per_L = (hmax-h0)/100;
        
    case 105
        % Based on Set#104
        expnt = 2;
        h0 = 320;
        hmax = 172;
        L_off = -20;
        L0   = 5;
        Lmax = 85;
        maxc = []; %[];
        
        h_per_L = (hmax-h0)/100;
        
    otherwise
        error('Unfamiliar method');
        
end

L = L0:Lmax;
Lmid = (L0+Lmax)/2;

h = h0 + h_per_L * ((L-L0) + L_off - (L-Lmid).^2 * L_off/(Lmax-Lmid)^2);

% len = sqrt((L(end)-L(1))^2 + h_per_L*(L(end)-L(1))^2);
% h = (L_off - (L-Lmid).^2 * L_off/(100-Lmid)^2);
% rotate (and scale) so that [1 0] -> [1 h_per_L]
% which means that [0 1] -> [-h_per_L 1]
% theta = atan(h_per_L);
% sf = 1/cos(theta);


h1 = round(h);
c1 = nan(size(L));
for i=1:length(L)
    tmp = g.lch(g.lch(:,1)==L(i),:); % Find gamut points at this L value
    [ignore,I] = min(abs(tmp(:,3)-h1(i))); % Find the gamut point closest in hue
    c1(i) = tmp(I,2);
end
maxc1 = min(c1(round(length(L)/2)+(-5:5))); % Min chroma of middle 10 values

c = 1 - abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(Lmax-Lmid,Lmid-L0))).^expnt) / abs(Lmid.^expnt);
if ~isempty(maxc)
    c = c*(maxc-c0);
end

% Go from cylindrical Lch co-ords to cartesian Lab
a = c.*cos(h/360*(2*pi));
b = c.*sin(h/360*(2*pi));

% Turn Lab into sRGB values
Lab = [L' a' b'];
rgb = gd_lab2rgb(Lab, use_uplab);

% Check for points out of gamut
[TF,P2] = isingamut(Lab,g,'Lab');

P2C = sqrt(P2(:,2).^2+P2(:,3).^2);
maxc2 = min(P2C(round(length(L)/2)+(-5:5)));
maxc3 = min(P2C./c');
if ~isempty(maxc)
    maxc3 = maxc3*maxc;
end

fprintf('Method #%d maxc1 = %.2f\tmaxc2 = %.2f\tmaxc3 = %.2f\n',method,maxc1,maxc2,maxc3);

if ~isempty(maxc)
    fprintf('Method #%d has %d of %d in gamut.\n',method,sum(TF),length(TF));
else
    % Redo with maxc3
    disp('Redoing with maxc3');
    c = (maxc3-c0) * (1 - abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(Lmax-Lmid,Lmid-L0))).^expnt) / abs(Lmid.^expnt));

    % Go from cylindrical Lch co-ords to cartesian Lab
    a = c.*cos(h/360*(2*pi));
    b = c.*sin(h/360*(2*pi));

    % Turn Lab into sRGB values
    Lab = [L' a' b'];
    rgb = gd_lab2rgb(Lab, use_uplab);

    % Check for points out of gamut
    [TF,P2] = isingamut(Lab,g,'Lab');
    fprintf('Method #%d now has %d of %d in gamut.\n',method,sum(TF),length(TF));

    P2C = sqrt(P2(:,2).^2+P2(:,3).^2);
end



% First and last colors may be out due to floating point error
% others should not be
% Lab(~TF,:) = P2(~TF,:);

% Plot map and its derivation
if dbg
img = repmat(rgb,[1 1 20]);
img = permute(img,[1 3 2]);

figure;
imagesc(img);
axis xy;
title(sprintf('expnt = %.3f; h0 = %.3f; h/L = %.3f; L.off = %.3f',...
    expnt, h0, h_per_L, L_off));

figure;
hold on;
plot(c1, L, 'Color', [0 .8 0]);
plot(c , L, 'k');
plot(sqrt(P2(:,2).^2+P2(:,3).^2), P2(:,1), 'r');
title(sprintf('expnt = %.3f; h0 = %.3f; h/L = %.3f; L.off = %.3f',...
    expnt, h0, h_per_L, L_off));
xlabel('chroma');
ylabel('Lightness');

figure;
hold on;
plot(L,c,'k');
plot(L,h,'r');
title(sprintf('expnt = %.3f; h0 = %.3f; h/L = %.3f; L.off = %.3f',...
    expnt, h0, h_per_L, L_off));
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
end