% Find seperation distance for a set of different hues
% (Blue-white-red style)
clear;
% close all;

%% Parameters

% CIELab bwr: 296, 40
% UPLab bwr:  309, 54.75
% handpicked_hue1 = 309; %[]; % 290; % Blue
% handpicked_hue2 = 54.75; %[]; % 41;  % Red

% hue1_range = 307:.25:312; %200:320; %260:315; %0:359; %260:315;
% hue2_range =  52:.25:56; % 20:90;  % 10:70;  %0:359; % 10:70 ;

% % bwr search
% hue1_range = 290:1:307;
% hue2_range =  39:1:60;

% hue1_range = 296;
% hue2_range =  40;

% hue1_range = 303;
% hue2_range =  41;

hue1_range = 296;
hue2_range =  41;
Lmid  = 91;
Ledg  = 15;

% Cool
% hue1_range = 103;
% hue2_range = 196;
% Lmid  = 15;
% Ledg  = 90;

% Curve parameters
use_uplab = false;
npoints = 100;     % number of points to use in test series
typ   = 'sin';     % 'pow' or 'sin'
c0rel = 0;
nLmax = 401;

switch typ
    case 'sin'
        expnt = 1;
    case 'pow'
        expnt = 2;
end

%% Main body

g = fetch_cielchab_gamut('srgb', [], [], use_uplab);


all_maxc = nan( length(hue1_range), length(hue2_range) );

li_L   = g.lchmesh.Lvec>=min(Ledg,Lmid) & g.lchmesh.Lvec<=max(Ledg,Lmid);
jointL = g.lchmesh.Lvec(li_L)';
L = jointL;

% Generate a set of higher precision L values as candidates for
% maximum chroma point
% This has to be in the lower half of the curve because C is tighter in
% the middle than the edge. Otherwise it wouldn't be called the middle.
LmaxCs = linspace((Lmid+Ledg)/2, Ledg, nLmax);

% Go through all hues
for ih1=1:length(hue1_range)
    
    hue1 = hue1_range(ih1);
    gh1  = [jointL g.lchmesh.cgrid(g.lchmesh.hvec==hue1,li_L)' g.lchmesh.hgrid(g.lchmesh.hvec==hue1,li_L)'];
    
    for ih2=1:length(hue2_range)
        
        hue2 = hue2_range(ih2);
        gh2  = [jointL g.lchmesh.cgrid(g.lchmesh.hvec==hue2,li_L)' g.lchmesh.hgrid(g.lchmesh.hvec==hue2,li_L)'];
        
        maxC = min(gh1(:,2),gh2(:,2));
        
        switch typ
            case 'sin'
                % Try a cosine curve from Ledg to Lmid with a start at Lmid and
                % a peak (pi/2) at LmaxCs(i)
                c = c0rel + (1-c0rel) * cos(pi* bsxfun(@rdivide,bsxfun(@minus,L,LmaxCs),2*abs(Lmid-LmaxCs)) ).^expnt;
            case 'pow'
                c = 1 - (1-c0rel) * bsxfun(@rdivide,bsxfun(@minus,L,LmaxCs),abs(Lmid-LmaxCs)).^expnt;
                c = max(0,c);
            otherwise
                error('Unfamiliar type');
        end
        % Divide gamut edge by generated curve C values to find maxc,
        % which is our scaling factor. Min to get best obtainable for
        % individual candidate. Max to get best candidate maxc.
        my_maxc = max(min(bsxfun(@rdivide,maxC,c)));
        
        all_maxc(ih1,ih2) = my_maxc;
        
    end
end


%%

figure;
imagesc(hue2_range, hue1_range, all_maxc);
colormap(clab_hot);
colorbar;
ylabel('Start hue');
xlabel('End hue');
title('Maximum chroma');


%%

params = struct;

params.use_uplab = use_uplab;
params.n     = 32;
params.typ   = typ;
params.expnt = expnt;
params.Lmid  = Lmid;
params.Ledg  = Ledg;


[M,I] = max(all_maxc(:));
casediscr = sprintf('Max Chroma = %.2f',M);

[mih1,mih2] = ind2sub(size(all_maxc),I);
hue1 = hue1_range(mih1);
hue2 = hue2_range(mih2);


% Need to work out best curve shape for these parameters
gh1  = [jointL g.lchmesh.cgrid(g.lchmesh.hvec==hue1,li_L)' g.lchmesh.hgrid(g.lchmesh.hvec==hue1,li_L)'];
gh2  = [jointL g.lchmesh.cgrid(g.lchmesh.hvec==hue2,li_L)' g.lchmesh.hgrid(g.lchmesh.hvec==hue2,li_L)'];
maxC = min(gh1(:,2),gh2(:,2));

% Increase precision on Lmaxc even further
LmaxCs = linspace((params.Lmid+params.Ledg)/2, params.Ledg, 4*nLmax);
switch typ
    case 'sin'
        % Try a cosine curve from Ledg to Lmid with a start at Lmid and
        % a peak (pi/2) at LmaxCs(i)
        c = c0rel + (1-c0rel) * cos(pi* bsxfun(@rdivide,bsxfun(@minus,L,LmaxCs),2*abs(Lmid-LmaxCs)) ).^expnt;
    case 'pow'
        c = 1 - (1-c0rel) * bsxfun(@rdivide,bsxfun(@minus,L,LmaxCs),abs(Lmid-LmaxCs)).^expnt;
        c = max(0,c);
    otherwise
        error('Unfamiliar type');
end
% Divide gamut edge by generated curve C values to find maxc,
% which is our scaling factor. Min to get best obtainable for
% individual candidate. Max to get best candidate maxc.
[my_maxc,I] = max(min(bsxfun(@rdivide,maxC,c)));
my_Lmaxc = LmaxCs(I);


params.h1    = hue1;
params.h2    = hue2;
params.Lmaxc = my_Lmaxc;
params.maxc  = my_maxc;
params.c0    = c0rel*params.maxc;

rgb = makecmap_AwpBcurve(params, true);

