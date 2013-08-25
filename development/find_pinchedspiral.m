
% handpicked_hue1 = 13;
% handpicked_hue2 = 87;
% 
% handpicked_hue1 = 21;
% handpicked_hue2 = 75;


% All up to one complete twist
% hue1_range = 0:359;
% hue2_range = -360:720;

% % All up to one complete twist, reduced sampling
% hue1_range =    0:2:359;
% hue2_range = -360:2:720;
% 
% % Hot seach space
% hue1_range = -45:1:45;
% hue2_range =  70:1:130;

% Hot seach space (tighter)
hue1_range =  0:1:30;
hue2_range = 70:1:110;



% Curve parameters
use_uplab = false;
npoints = 100;  % number of points to use in test series
expnt = 2.3;    % exponent to use for chroma curve
L_off = 0;      % offset to exponent (not for sine curve)
c0 = 0.3;       % proportion of maxc to add to all chroma (70<maxc<80)
typ = 'pow';    % 'pow' or 'sin'
Lmin = 5;       % minimum lightness (start)
Lmax = 95;      % maximum lightness (end)

% Computationally defined parameters
L = linspace(Lmin,Lmax,npoints);
Lmid = (Lmin+Lmax)/2;


g = fetch_cielchab_gamut('srgb', [], [], use_uplab);

all_maxc = nan(length(hue1_range),length(hue2_range));
all_eucl = nan(length(hue1_range),length(hue2_range));

for ih1=1:length(hue1_range)
    for ih2=1:length(hue2_range);
        h1 = hue1_range(ih1);
        h2 = hue2_range(ih2);
        h_per_L = (h2-h1)/(Lmax-Lmin);

        h = h1 + h_per_L * ((L-Lmin) + L_off - (L-Lmid).^2 * L_off/(Lmax-Lmid)^2);
        h = mod(h,360);
        
        switch typ
            case 'sin'
                c = c0 + (1-c0) * sin(pi* (L-Lmin)/(Lmax-Lmin) ).^expnt;
            case 'pow'
                c = 1 - (1-c0) * abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(Lmax-Lmid,Lmid-Lmin))).^expnt) / abs(Lmid.^expnt);
            otherwise
                error('Unfamiliar type');
        end

        Lch = [L' c' h'];

        % Check for points out of gamut
        [TF,P2] = isingamut(Lch,g,'Lch');
        
        P2C = P2(:,2);
        maxc = min(P2C./c');
        all_maxc(ih1,ih2) = maxc;
        
        % Euclidean length
        c = c * maxc;
        a = c.*cosd(h);
        b = c.*sind(h);
        Lab = [L' a' b'];
        Lab_dif = diff(Lab,1,1);
        Lab_sep = sqrt(sum(Lab_dif.^2,2));
        Lab_len = sum(Lab_sep,1);
        
        all_eucl(ih1,ih2) = Lab_len;
    end
end

all_both = all_maxc/max(all_maxc(:))+all_eucl/max(all_eucl(:));

%%

figure;
imagesc(hue2_range, hue1_range, all_maxc);
colormap(cie_hot);
colorbar;
ylabel(sprintf('Start hue at L=%.1f',Lmin));
xlabel(sprintf('End hue at L=%.1f',Lmax));
title('Maximum chroma');

figure;
imagesc(hue2_range, hue1_range, all_eucl);
colormap(cie_hot);
colorbar;
ylabel(sprintf('Start hue at L=%.1f',Lmin));
xlabel(sprintf('End hue at L=%.1f',Lmax));
title('Euclidean curve length');

figure;
imagesc(hue2_range,hue1_range,all_both);
colormap(cie_hot);
colorbar;
ylabel(sprintf('Start hue at L=%.1f',Lmin));
xlabel(sprintf('End hue at L=%.1f',Lmax));
title('Sum');

%%
% return;

%%

% h1 = 275;
% h2 = 360;

% h1 = 290;
% h2 = 380;



for i=1:3
    switch i
        case 0
            if ~exist('handpicked_hue1','var') || ~exist('handpicked_hue2','var')
                continue;
            end
            h1 = handpicked_hue1;
            h2 = handpicked_hue2;
            casediscr = 'Handpicked';
        case 1
            [mm,mI] = max(all_maxc(:));
            [mih1,mih2] = ind2sub(size(all_maxc),mI);
            h1 = hue1_range(mih1);
            h2 = hue2_range(mih2);
            casediscr = sprintf('Max Chroma = %.2f',mm);
        case 2
            [mm,mI] = max(all_eucl(:));
            [mih1,mih2] = ind2sub(size(all_maxc),mI);
            h1 = hue1_range(mih1);
            h2 = hue2_range(mih2);
            casediscr = sprintf('Max Euclidean Length = %.2f',mm);
        case 3
            [mm,mI] = max(all_both(:));
            [mih1,mih2] = ind2sub(size(all_maxc),mI);
            h1 = hue1_range(mih1);
            h2 = hue2_range(mih2);
            casediscr = sprintf('Max weighted sum = %.2f',mm);
    end


h_per_L = (h2-h1)/(Lmax-Lmin);

h = h1 + h_per_L * ((L-Lmin) + L_off - (L-Lmid).^2 * L_off/(Lmax-Lmid)^2);
h = mod(h,360);

switch typ
    case 'sin'
        c = c0 + (1-c0) * sin(pi* (L-Lmin)/(Lmax-Lmin) ).^expnt;
    case 'pow'
        c = 1 - (1-c0) * abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(Lmax-Lmid,Lmid-Lmin))).^expnt) / abs(Lmid.^expnt);
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
title(sprintf('%.2f %.2f: %s',h1,h2,casediscr));

plot_labcurve_rgbgamut(Lab, use_uplab);
title(sprintf('%.2f %.2f: %s',h1,h2,casediscr));

end

%%

% h1 = handpicked_hue1;
% h2 = handpicked_hue2;

return;
%%
[m1,I1] = max(all_both,[],2);

for ih1=1:3:length(hue1_range)
    ih2 = I1(ih1);
    h1 = hue1_range(ih1);
    h2 = hue2_range(ih2);
    casediscr = '';

h_per_L = (h2-h1)/(Lmax-Lmin);

h = h1 + h_per_L * ((L-Lmin) + L_off - (L-Lmid).^2 * L_off/(Lmax-Lmid)^2);
h = mod(h,360);

switch typ
    case 'sin'
        c = c0 + (1-c0) * sin(pi* (L-Lmin)/(Lmax-Lmin) ).^expnt;
    case 'pow'
        c = 1 - (1-c0) * abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(Lmax-Lmid,Lmid-Lmin))).^expnt) / abs(Lmid.^expnt);
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
title(sprintf('%.2f %.2f: %s',h1,h2,casediscr));

% plot_labcurve_rgbgamut(Lab, use_uplab);
% title(sprintf('%.2f %.2f: %s',h1,h2,casediscr));

end