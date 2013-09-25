% Find seperation distance for a set of different hues
% (Blue-white-red style)
clear all;
% close all;

%% Parameters

use_uplab = false;

% CIELab bwr: 296, 40
% UPLab bwr:  309, 54.75
handpicked_hue1 = 309; %[]; % 290; % Blue
handpicked_hue2 = 54.75; %[]; % 41;  % Red

hue1_range = 307:.25:312; %200:320; %260:315; %0:359; %260:315;
hue2_range =  52:.25:56; % 20:90;  % 10:70;  %0:359; % 10:70 ;

% hue1_range =  40;
% hue2_range = 296;

% bwr
hue1_range = 285:307;
hue2_range =  30:50;
via_black = false;

% bwr tight
hue1_range = 291:.5:306;
hue2_range =  33:.5:50;
via_black = false;


% Curve parameters
use_uplab = false;
npoints = 100;     % number of points to use in test series
typ = 'sin';       % 'pow' or 'sin'
c0 = 0;
% Lmin     =  0; % curve shape
Lmax     = 94; % curve shape
spotLmin = 15; % displayed colour limit
spotLmax = Lmax; % displayed colour limit
ncurve   = 401;

switch typ
    case 'sin'
        expnt = 1;
    case 'pow'
        expnt = 2;
end

%% Main body

g = fetch_cielchab_gamut('srgb', [], [], use_uplab);


all_maxc = nan( length(hue1_range), length(hue2_range), ncurve );

% Computationally defined parameters
% L = linspace(Lmin,Lmax,npoints);

if via_black
    % Sets Lmax
    Lcurve = linspace(spotLmax, spotLmax-Lmin, ncurve);
else
    % Sets Lmin
    Lcurve = linspace(spotLmin, spotLmin-Lmax, ncurve);
end

li_L   = g.lchmesh.Lvec>=spotLmin & g.lchmesh.Lvec<=spotLmax;
jointL = g.lchmesh.Lvec(li_L)';
L = jointL;

% Go through all hue
for ih1=1:length(hue1_range)
    
    hue1 = hue1_range(ih1);
    li_h = g.lchmesh.hvec==hue1;
    gh1  = [jointL g.lchmesh.cgrid(li_h,li_L)' g.lchmesh.hgrid(li_h,li_L)'];
    
    for ih2=1:length(hue2_range)
        
        hue2 = hue2_range(ih2);
        li_h = g.lchmesh.hvec==hue2;
        gh2  = [jointL g.lchmesh.cgrid(li_h,li_L)' g.lchmesh.hgrid(li_h,li_L)'];
        
        jointC = min(gh1(:,2),gh2(:,2));
        
        % Try different curve steepness to see how well we can do
        for icurve=1:ncurve
            
            if via_black
                Lmax = Lcurve(icurve);
            else
                Lmin = Lcurve(icurve);
            end
            Lmid = (Lmin+Lmax)/2;
            
            switch typ
                case 'sin'
                    c = c0 + (1-c0) * sin(pi* (L-Lmin)/(Lmax-Lmin) ).^expnt;
                case 'pow'
                    c = 1 - (1-c0) * abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(Lmax-Lmid,Lmid-Lmin))).^expnt) / abs(Lmid.^expnt);
                    c = max(0,c);
                otherwise
                    error('Unfamiliar type');
            end
            
            % Check for points out of gamut
            maxc = min(jointC./c);
            all_maxc(ih1,ih2,icurve) = maxc;
            
        end
        
    end
end

[all_maxc,Icurvemax] = max(all_maxc,[],3);


%%

figure;
imagesc(hue2_range, hue1_range, all_maxc);
colormap(clab_hot);
colorbar;
ylabel('Start hue');
xlabel('End hue');
title('Maximum chroma');

%%

[mm,mIh] = max(all_maxc(:));

if via_black
    Lmax = Lcurve(Icurvemax(mIh));
else
    Lmin = Lcurve(Icurvemax(mIh));
end
Lmid = (Lmin+Lmax)/2;

casediscr = sprintf('Max Chroma = %.2f',mm);
[mih1,mih2] = ind2sub(size(all_maxc),mIh);

h1 = hue1_range(mih1);
h2 = hue2_range(mih2);

L = linspace(spotLmin,spotLmax,floor(npoints/2)+1);

switch typ
    case 'sin'
        c = c0 + (1-c0) * sin(pi* (L-Lmin)/(Lmax-Lmin) ).^expnt;
    case 'pow'
        c = 1 - (1-c0) * abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(Lmax-Lmid,Lmid-Lmin))).^expnt) / abs(Lmid.^expnt);
        c = max(0,c);
    otherwise
        error('Unfamiliar type');
end

% Check for points out of gamut
maxc = all_maxc(mih1,mih2);

c = c*maxc;

Lch1 = [L' c' repmat(h1,size(L))'];
Lch2 = [flipud(L') flipud(c') repmat(h2,size(L))'];

if c0==0
    Lch1 = Lch1(1:end-1,:);
end
lch = [Lch1;Lch2];

lab = [lch(:,1) lch(:,2).*cosd(lch(:,3)) lch(:,2).*sind(lch(:,3))];
% lch = [lab(1) sqrt(lab(2)^2+lab(3)^2) mod(atan2(lab(3),lab(2))/pi*180,360)];

rgb = hard_lab2rgb(lab, use_uplab);


img = repmat(rgb,[1 1 20]);
img = permute(img,[1 3 2]);
figure;
imagesc(img);
% title(sprintf('%.2f %.2f; %.2f %.2f; %.2f: %s',h1,h2,Lmin,Lmax,expnt,casediscr));

plot_labcurve_rgbgamut(lab, use_uplab);
% title(sprintf('%.2f %.2f; %.2f %.2f; %.2f: %s',h1,h2,Lmin,Lmax,expnt,casediscr));


% Debug
jointL = g.lchmesh.Lvec(li_L)';
li_h = g.lchmesh.hvec==h1;
gh1  = [jointL g.lchmesh.cgrid(li_h,li_L)' g.lchmesh.hgrid(li_h,li_L)'];
li_h = g.lchmesh.hvec==h2;
gh2  = [jointL g.lchmesh.cgrid(li_h,li_L)' g.lchmesh.hgrid(li_h,li_L)'];
jointC = min(gh1(:,2),gh2(:,2));

figure; set(gca,'Color',[.467 .467 .467]); hold on; box on;
plot(gh1(:,2),gh1(:,1),'b-');
plot(gh2(:,2),gh2(:,1),'r-');
plot(c,L,'ks-');


%%
return;
% Should pick points which are equispaced instead
%%
plot_labcurve_rgbgamut(hard_rgb2lab(cbrewer('div','RdBu',11)));

