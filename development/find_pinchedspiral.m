close all;

%%

% handpicked_hue1 = 13;
% handpicked_hue2 = 87;
% 
% handpicked_hue1 = 21;
% handpicked_hue2 = 75;

% All max one twist maximum euclidian length
handpicked_hue1 = 252;
handpicked_hue2 = 564;

% All max one twist maximum chroma


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

% % Hot seach space (tighter)
% hue1_range =  0:1:30;
% hue2_range = 70:1:110;

% % Hot seach space (tighter)
% hue1_range =  0:1:30;
% hue2_range = 70:1:106;

% % Hot seach space (tighter still)
% hue1_range =  7:1:25;
% hue2_range = 95:1:102;
% % Outcome was h1=9;h2=95;Lmin=4;Lmax=96;expnt=2;

% Cold search space
hue1_range = 300:2:360;
hue2_range = 210:2:260;

% % Purple
% hue1_range = 290:2:330;
% hue2_range = 310:2:350;
% % Outcome h1=308.5;h2=331;Lmin=4;Lmax=92;expnt=2;

hue1_range = 301:.5:312;
hue2_range = 328:.5:335;

% hue1_range = 300:2:360;
% hue2_range = 200:2:300;


% Curve parameters
use_uplab = false;
npoints = 100;     % number of points to use in test series
expnt_range = 2;   % exponents to use for chroma curve
L_off = 0;         % offset to exponent (not for sine curve)
c0 = 0.025;         % proportion of maxc to add to all chroma (maxc\approx75, 67<maxc<77)
typ = 'pow';       % 'pow' or 'sin'
Lmin_range =  4;   % minimum lightness (start)
Lmax_range = 92;   % maximum lightness (end)



g = fetch_cielchab_gamut('srgb', [], [], use_uplab);

all_maxc = nan( length(hue1_range), length(hue2_range), length(Lmin_range), length(Lmax_range), length(expnt_range));
all_eucl = nan( length(hue1_range), length(hue2_range), length(Lmin_range), length(Lmax_range), length(expnt_range));

for iXpt=1:length(expnt_range)
    expnt = expnt_range(iXpt);
    
    for iLmin=1:length(Lmin_range)
        Lmin = Lmin_range(iLmin);
    
        for iLmax=1:length(Lmax_range)
            Lmax = Lmax_range(iLmax);
        
            % Computationally defined parameters
            L = linspace(Lmin,Lmax,npoints);
            Lmid = (Lmin+Lmax)/2;

            for ih1=1:length(hue1_range)
                h1 = hue1_range(ih1);

                for ih2=1:length(hue2_range);
                    h2 = hue2_range(ih2);
                    h_per_L = (h2-h1)/(Lmax-Lmin);

                    h = h1 + h_per_L * ((L-Lmin) + L_off - (L-Lmid).^2 * L_off/(Lmax-Lmid)^2);
                    h = mod(h,360);

                    switch typ
                        case 'sin'
                            c = c0 + (1-c0) * sin(pi* (L-Lmin)/(Lmax-Lmin) ).^expnt;
                        case 'pow'
                            c = 1 - (1-c0) * abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(Lmax-Lmid,Lmid-Lmin))).^expnt) / abs(Lmid.^expnt);
                            c = max(0,c);
                        otherwise
                            error('Unfamiliar type');
                    end

                    Lch = [L' c' h'];

                    % Check for points out of gamut
                    [TF,P2] = isingamut(Lch,g,'Lch');

                    P2C = P2(:,2);
                    maxc = min(P2C./c');
                    all_maxc(ih1,ih2,iLmin,iLmax,iXpt) = maxc;

                    % Euclidean length
                    c = c * maxc;
                    a = c.*cosd(h);
                    b = c.*sind(h);
                    Lab = [L' a' b'];
                    Lab_dif = diff(Lab,1,1);
                    Lab_sep = sqrt(sum(Lab_dif.^2,2));
                    Lab_len = sum(Lab_sep,1);

                    all_eucl(ih1,ih2,iLmin,iLmax,iXpt) = Lab_len;
                end
            end
        end
    end
end

all_both = 0.7*all_maxc/max(all_maxc(:)) + 0.3*all_eucl/max(all_eucl(:));

all_maxc2 = max(max(max(all_maxc,[],5),[],4),[],3);
all_eucl2 = max(max(max(all_eucl,[],5),[],4),[],3);
all_both2 = max(max(max(all_both,[],5),[],4),[],3);

%%

figure;
imagesc(hue2_range, hue1_range, all_maxc2);
colormap(clab_hot);
colorbar;
ylabel('Start hue');
xlabel('End hue');
title('Maximum chroma');

figure;
imagesc(hue2_range, hue1_range, all_eucl2);
colormap(clab_hot);
colorbar;
ylabel('Start hue');
xlabel('End hue');
title('Euclidean curve length');

figure;
imagesc(hue2_range,hue1_range,all_both2);
colormap(clab_hot);
colorbar;
ylabel('Start hue');
xlabel('End hue');
title('Sum');

%%
if length(Lmin_range)>1 || length(Lmax_range)>1

all_maxc2L = max(max(max(all_maxc,[],5),[],1),[],2);
all_eucl2L = max(max(max(all_eucl,[],5),[],1),[],2);

figure;
imagesc(Lmin_range, Lmax_range,  squeeze(all_maxc2L)');
colormap(clab_hot);
colorbar;
xlabel('Lmin');
ylabel('Lmax');
title('Maximum chroma');

figure;
imagesc(Lmin_range, Lmax_range,  squeeze(all_eucl2L)');
colormap(clab_hot);
colorbar;
xlabel('Lmin');
ylabel('Lmax');
title('Euclidean curve length');

end

%%
if length(expnt_range)>1

all_maxc2X = max(max(max(max(all_maxc,[],4),[],3),[],2),[],1);
all_eucl2X = max(max(max(max(all_eucl,[],4),[],3),[],2),[],1);

figure;
plot(expnt_range, squeeze(all_maxc2X));
xlabel('Exponent');
ylabel('Maximum chroma');

figure;
plot(expnt_range, squeeze(all_eucl2X));
xlabel('Exponent');
ylabel('Euclidean curve length');

end

%%

% h1 = 275;
% h2 = 360;

% h1 = 290;
% h2 = 380;



for i=1:3
    switch i
%         case 0
%             if ~exist('handpicked_hue1','var') || ~exist('handpicked_hue2','var')
%                 continue;
%             end
%             h1 = handpicked_hue1;
%             h2 = handpicked_hue2;
%             casediscr = 'Handpicked';
        case 1
%             [mm,mIh] = max(all_maxc2(:));
%             casediscr = sprintf('Max Chroma = %.2f',mm);
%             [mih1,mih2] = ind2sub(size(all_maxc2),mIh);
%             mIL = find(all_maxc(mih1,mih2,:,:)>=mm);
%             [~,~,miLmin,miLmax] = ind2sub(size(all_maxc(mih1,mih2,:,:)),mIL);
            [mm,mIh] = max(all_maxc(:));
            casediscr = sprintf('Max Chroma = %.2f',mm);
            [mih1,mih2,miLmin,miLmax,miXpt] = ind2sub(size(all_maxc),mIh);
            
        case 2
%             [mm,mIh] = max(all_eucl2(:));
%             casediscr = sprintf('Max Euclidean Length = %.2f',mm);
%             [mih1,mih2] = ind2sub(size(all_eucl2),mIh);
%             mIL = find(all_eucl2(mih1,mih2,:,:)>=mm);
%             [~,~,miLmin,miLmax] = ind2sub(size(all_eucl2(mih1,mih2,:,:)),mIL);
            [mm,mIh] = max(all_eucl(:));
            casediscr = sprintf('Max Euclidean Length = %.2f',mm);
            [mih1,mih2,miLmin,miLmax,miXpt] = ind2sub(size(all_eucl),mIh);
        case 3
%             [mm,mIh] = max(all_both2(:));
%             casediscr = sprintf('Max weighted sum = %.2f',mm);
%             [mih1,mih2] = ind2sub(size(all_both2),mIh);
%             mIL = find(all_both2(mih1,mih2,:,:)>=mm);
%             [~,~,miLmin,miLmax] = ind2sub(size(all_both2(mih1,mih2,:,:)),mIL);
            [mm,mIh] = max(all_both(:));
            casediscr = sprintf('Max weighted sum = %.2f',mm);
            [mih1,mih2,miLmin,miLmax,miXpt] = ind2sub(size(all_both),mIh);
    end
    
    h1 = hue1_range(mih1);
    h2 = hue2_range(mih2);
    Lmin = Lmin_range(miLmin);
    Lmax = Lmax_range(miLmax);
    expnt = expnt_range(miXpt);
    
    
    L = linspace(Lmin,Lmax,npoints);
    Lmid = (Lmin+Lmax)/2;

    h_per_L = (h2-h1)/(Lmax-Lmin);

    h = h1 + h_per_L * ((L-Lmin) + L_off - (L-Lmid).^2 * L_off/(Lmax-Lmid)^2);
    h = mod(h,360);

    switch typ
        case 'sin'
            c = c0 + (1-c0) * sin(pi* (L-Lmin)/(Lmax-Lmin) ).^expnt;
        case 'pow'
            c = 1 - (1-c0) * abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(Lmax-Lmid,Lmid-Lmin))).^expnt) / abs(Lmid.^expnt);
            c = max(0,c);
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
    title(sprintf('%.2f %.2f; %.2f %.2f; %.2f: %s',h1,h2,Lmin,Lmax,expnt,casediscr));

    plot_labcurve_rgbgamut(Lab, use_uplab);
    title(sprintf('%.2f %.2f; %.2f %.2f; %.2f: %s',h1,h2,Lmin,Lmax,expnt,casediscr));

end

%%

% h1 = handpicked_hue1;
% h2 = handpicked_hue2;
% Lmin = Lmin_range(1);
% Lmax = Lmax_range(1);
% expnt = expnt_range(1);

%%

for iXpt=1:length(expnt_range)
    expnt = expnt_range(iXpt);
    
    ix_all_both = all_both(:,:,:,:,iXpt);
    [mm,mIh] = max(ix_all_both(:));
    casediscr = sprintf('Max weighted sum = %.2f',mm);
    [mih1,mih2,miLmin,miLmax,miXpt] = ind2sub(size(ix_all_both),mIh);

    h1 = hue1_range(mih1);
    h2 = hue2_range(mih2);
    Lmin = Lmin_range(miLmin);
    Lmax = Lmax_range(miLmax);
    
L = linspace(Lmin,Lmax,npoints);
Lmid = (Lmin+Lmax)/2;

h_per_L = (h2-h1)/(Lmax-Lmin);

h = h1 + h_per_L * ((L-Lmin) + L_off - (L-Lmid).^2 * L_off/(Lmax-Lmid)^2);
h = mod(h,360);

switch typ
    case 'sin'
        c = c0 + (1-c0) * sin(pi* (L-Lmin)/(Lmax-Lmin) ).^expnt;
    case 'pow'
        c = 1 - (1-c0) * abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(Lmax-Lmid,Lmid-Lmin))).^expnt) / abs(Lmid.^expnt);
        c = max(0,c);
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

%%
figure(round(100*expnt));
imagesc(hue2_range,hue1_range,all_both2);
colormap(rgb);
colorbar;
ylabel('Start hue');
xlabel('End hue');
title(sprintf('%.2f %.2f; %.2f %.2f; %.2f: Sum (%.2f, %.2f)',...
    h1,h2,Lmin,Lmax,expnt,all_maxc(mih1,mih2,miLmin,miLmax,iXpt),all_eucl(mih1,mih2,miLmin,miLmax,iXpt)));

% figure(101); hold on;
% plot(L,h);
% xlabel('L'); ylabel('hue');
% 
% figure(102); hold on;
% plot(L,c);
% xlabel('L'); ylabel('chroma');
% 
% Lab_dif = sqrt(sum(diff(Lab,1,1).^2,2));
% figure(103); hold on;
% plot(L(2:end),Lab_dif);
% xlabel('L'); ylabel('Euclidean point sep');

end

%%
return;
%%
[m1,I1] = max(all_both2,[],2);

for ih1=1:3:length(hue1_range)
    ih2 = I1(ih1);
    h1 = hue1_range(ih1);
    h2 = hue2_range(ih2);
    casediscr = '';

L = linspace(Lmin,Lmax,npoints);
Lmid = (Lmin+Lmax)/2;

h_per_L = (h2-h1)/(Lmax-Lmin);

h = h1 + h_per_L * ((L-Lmin) + L_off - (L-Lmid).^2 * L_off/(Lmax-Lmid)^2);
h = mod(h,360);

switch typ
    case 'sin'
        c = c0 + (1-c0) * sin(pi* (L-Lmin)/(Lmax-Lmin) ).^expnt;
    case 'pow'
        c = 1 - (1-c0) * abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(Lmax-Lmid,Lmid-Lmin))).^expnt) / abs(Lmid.^expnt);
        c = max(0,c);
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

%%
return;
%%

h1 = handpicked_hue1;
h2 = handpicked_hue2;
h1=308.5;h2=331;Lmin=4;Lmax=92;expnt=2;

% ih1 = find(h1==hue1_range,1);
% ih2 = find(h2==hue2_range,1);

    Lmin = Lmin_range(miLmin);
    Lmax = Lmax_range(miLmax);
    expnt = expnt_range(miXpt);
    
    
L = linspace(Lmin,Lmax,npoints);
Lmid = (Lmin+Lmax)/2;

h_per_L = (h2-h1)/(Lmax-Lmin);

h = h1 + h_per_L * ((L-Lmin) + L_off - (L-Lmid).^2 * L_off/(Lmax-Lmid)^2);
h = mod(h,360);

switch typ
    case 'sin'
        c = c0 + (1-c0) * sin(pi* (L-Lmin)/(Lmax-Lmin) ).^expnt;
    case 'pow'
        c = 1 - (1-c0) * abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(Lmax-Lmid,Lmid-Lmin))).^expnt) / abs(Lmid.^expnt);
        c = max(0,c);
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

% Euclidean length
Lab_dif = diff(Lab,1,1);
Lab_sep = sqrt(sum(Lab_dif.^2,2));
Lab_len = sum(Lab_sep,1);

rgb = hard_lab2rgb(Lab, use_uplab);

img = repmat(rgb,[1 1 20]);
img = permute(img,[1 3 2]);
figure;
imagesc(img);
title(sprintf('%.2f %.2f; %.2f %.2f; %.2f: Sum (%.2f, %.2f)',...
    h1,h2,Lmin,Lmax,expnt,maxc,Lab_len));

plot_labcurve_rgbgamut(Lab, use_uplab);
title(sprintf('%.2f %.2f; %.2f %.2f; %.2f: Sum (%.2f, %.2f)',...
    h1,h2,Lmin,Lmax,expnt,maxc,Lab_len));