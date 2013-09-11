% Find seperation distance for a set of different hues
% (Blue-white-red style)
clear all;
close all;

%% Parameters

use_uplab = false;

% CIELab bwr: 296, 40
% UPLab bwr:  309, 54.75
handpicked_hue1 = 309; %[]; % 290; % Blue
handpicked_hue2 = 54.75; %[]; % 41;  % Red

via_black = 0; % Go via near-white (0) or near-black (1)

hue1_range = 307:.25:312; %200:320; %260:315; %0:359; %260:315;
hue2_range =  52:.25:56; % 20:90;  % 10:70;  %0:359; % 10:70 ;

%% Main body

g = fetch_cielchab_gamut('srgb', [], [], use_uplab);


all_best_sep     = nan(length(hue1_range),length(hue2_range));
all_best_sep_m   = nan(length(hue1_range),length(hue2_range));
all_best_C_sep   = nan(length(hue1_range),length(hue2_range));
all_best_C_sep_m = nan(length(hue1_range),length(hue2_range));


for ih1=1:length(hue1_range)
    for ih2=1:length(hue2_range)
        hue1 = hue1_range(ih1);
        hue2 = hue2_range(ih2);

% % OLD METHOD: SLOW
% gh1 = g.lch(g.lch(:,3)==hue1,:);
% gh2 = g.lch(g.lch(:,3)==hue2,:);
% jointL = intersect(gh1(:,1),gh2(:,1));
% gh1 = gh1(ismember(gh1(:,1),jointL),:);
% gh2 = gh2(ismember(gh1(:,1),jointL),:);
% NEW METHOD
jointL = g.lchmesh.Lvec';
li = g.lchmesh.hvec==hue1;
gh1 = [jointL g.lchmesh.cgrid(li,:)' g.lchmesh.hgrid(li,:)'];
li = g.lchmesh.hvec==hue2;
gh2 = [jointL g.lchmesh.cgrid(li,:)' g.lchmesh.hgrid(li,:)'];


jointC = min(gh1(:,2),gh2(:,2));
[Cmax,I_Cmax] = max(jointC);
L_Cmax = gh2(I_Cmax,1);

switch via_black
    case 0
        valid_I = I_Cmax:length(jointC);
    case 1
        valid_I = 1:I_Cmax;
    otherwise
        error('Bad via colour');
end

% figure;hold on;set(gca,'Color',[.467 .467 .467]);box on;
% plot(gh1(:,1),gh1(:,2),'b-');
% plot(gh2(:,1),gh2(:,2),'r-');
% % plot(jointL,jointC,'kx');


m_segs = diff(jointC(valid_I))./diff(jointL(valid_I));
m_max = max(m_segs);
m_min = min(m_segs);

m_list = linspace(m_min,m_max,101);
initL = jointL(valid_I);

best_sep = 0;
best_sep_m = [];

best_C_sep = 0;
best_C_m = [];

% figure;hold on;set(gca,'Color',[.467 .467 .467]);box on;
% plot(jointL,jointC,'k-');

for i=1:length(m_list)
    m = m_list(i);
    % Given a gradient, compute the chroma values for each L
%     C = m*L+k
    % Since L=100 has C=0, curve is offset by
%     k = C-m*L
    if via_black
        initk = 0;
    else
        initk = -m*100;
    end
    initC = m*initL+initk;
    
    % Find min chroma offset needed so all datapoints are in gamut
    k = initk - max(initC - jointC(valid_I));
    
    % Find the intersection of this moved line with C=0
%     C = m*L+k
%     0 = m*Lend+k
%     Lend = -k/m
    Lend = -k/m;
    
    % Now we know the L start/end value needed to make this gradient work
    
    % Find the furthest we can draw the line before intersecting chroma
    % boundary
    if via_black
        moveL = jointL(jointL>=Lend);
    else
        moveL = jointL(jointL<=Lend);
    end
    moveC = m*moveL+k;
    if via_black
        I_start = find(moveC>jointC(jointL>=Lend),1,'first')-1;
        if I_start<0; continue; end;
    else
        I_start = find(moveC>jointC(jointL<=Lend),1,'last')+1;
        if I_start>length(moveC); continue; end;
    end
    Lstart = jointL(I_start);
    Cstart = moveC(I_start);
    
% %     figure(2);hold on;set(gca,'Color',[.467 .467 .467]);box on;
% %     plot(jointL,jointC,'k-');
%     plot(initL,initC,'r-');
%     plot(Lend,0,'ro');
%     plot(moveL,moveC,'b-');
%     plot(Lstart,Cstart,'bo');
    
    % Define metric of how good this gradient is
    % Euclidian distance from L start value and C max value
    sep_dist = sqrt((Lstart-Lend).^2+Cstart.^2);
    
    % If we have the best Euclidian distance, save this data
    if sep_dist>best_sep
        best_sep = sep_dist;
        best_sep_m = m;
    end
    
    % If we have the best Euclidian distance, save this data
    if Cstart>best_C_sep
        best_C_sep = Cstart;
        best_C_m = m;
    end
    
end

        all_best_sep(ih1,ih2)     = best_sep;
        all_best_sep_m(ih1,ih2)   = best_sep_m;
        all_best_C_sep(ih1,ih2)   = best_C_sep;
        all_best_C_sep_m(ih1,ih2) = best_C_m;
        
    end
end


%%

CMAP = cie_hot;

figure;
subplot(2,2,1);
imagesc(hue2_range,hue1_range,all_best_sep);
colormap(CMAP); colorbar;
% freezeColors; cbfreeze;
title('Best seperation distance');
subplot(2,2,2);
imagesc(hue2_range,hue1_range,all_best_sep_m);
colormap(CMAP); colorbar;
% freezeColors; cbfreeze;
title('gradient');
subplot(2,2,3);
imagesc(hue2_range,hue1_range,all_best_C_sep);
colormap(CMAP); colorbar;
% freezeColors; cbfreeze;
title('Best chroma distance');
subplot(2,2,4);
imagesc(hue2_range,hue1_range,all_best_C_sep_m);
colormap(CMAP); colorbar;
% freezeColors; cbfreeze;
title('gradient');

%% Best seperation

% m = best_sep_m;

[n_best_sep,I] = max(all_best_sep(:));
m = all_best_sep_m(I);
[ih1,ih2] = ind2sub(size(all_best_sep),I); %numel(all_best_sep));

hue1 = hue1_range(ih1);
hue2 = hue2_range(ih2);

gh1 = g.lch(g.lch(:,3)==hue1,:);
gh2 = g.lch(g.lch(:,3)==hue2,:);

jointL = intersect(gh1(:,1),gh2(:,1));
gh1 = gh1(ismember(gh1(:,1),jointL),:);
gh2 = gh2(ismember(gh1(:,1),jointL),:);

jointC = min(gh1(:,2),gh2(:,2));
[Cmax,I_Cmax] = max(jointC);
L_Cmax = gh2(I_Cmax,1);

switch via_black
    case 0
        valid_I = I_Cmax:length(jointC);
    case 1
        valid_I = 1:I_Cmax;
    otherwise
        error('Bad via colour');
end
initL = jointL(valid_I);

% % Plot the blue and red chroma curves to see how well they match
% figure;hold on;set(gca,'Color',[.467 .467 .467]);box on;
% plot(gh1(:,1),gh1(:,2),'b-');
% plot(gh2(:,1),gh2(:,2),'r-');
% % plot(jointL,jointC,'kx');
% title(sprintf('Best seperation distance (%d, %d) = %.2f',hue1,hue2,sep_dist));

    % Given a gradient, compute the chroma values for each L
%     C = m*L+k
    % Since L=100 has C=0, curve is offset by
%     k = C-m*L
    if via_black
        initk = 0;
    else
        initk = -m*100;
    end
    initC = m*initL+initk;
    
    % Find min chroma offset needed so all datapoints are in gamut
    k = initk - max(initC - jointC(valid_I));
    
    % Find the intersection of this moved line with C=0
%     C = m*L+k
%     0 = m*Lend+k
%     Lend = -k/m
    Lend = -k/m;
    
    % Now we know the L start value needed to make this gradient work
    
    % Find the furthest we can draw the line before intersecting chroma
    % boundary
    if via_black
        moveL = jointL(jointL>=Lend);
    else
        moveL = jointL(jointL<=Lend);
    end
    moveC = m*moveL+k;
    if via_black
        I_start = find(moveC>jointC(jointL>=Lend),1,'first')-1;
        if I_start<0; continue; end;
    else
        I_start = find(moveC>jointC(jointL<=Lend),1,'last')+1;
        if I_start>length(moveC); continue; end;
    end
    Lstart = jointL(I_start);
    Cstart = moveC(I_start);

    % Define metric of how good this gradient is
    % Euclidian distance from L start value and C max value
    sep_dist = sqrt((Lstart-Lend).^2+Cstart.^2);
    

% % Plot derivation of map
% figure;hold on;set(gca,'Color',[.467 .467 .467]);box on;
% plot(jointL,jointC,'k-');
% plot(initL,initC,'y-');
% plot(Lend,0,'yo');
% plot(moveL,moveC,'g-');
% plot(Lstart,Cstart,'go');
% title(sprintf('Chosen seperation distance (%d, %d) = %.2f',hue1,hue2,sep_dist));

% Plot the blue and red chroma curves to see how well they match
% WITH derivation of map
figure;hold on;set(gca,'Color',[.467 .467 .467]);box on;
plot(gh1(:,1),gh1(:,2),'b-');
plot(gh2(:,1),gh2(:,2),'r-');
plot(initL,initC,'y-');
plot(Lend,0,'yo');
plot(moveL,moveC,'g-');
plot(Lstart,Cstart,'go');
% plot(jointL,jointC,'kx');
title(sprintf('Chosen seperation distance (%d, %d) = %.2f',hue1,hue2,sep_dist));


% Plot sample colormap
neach = 32;
L1 = linspace(Lstart            , Lend, neach);
a1 = linspace(Cstart*cosd(hue1),    0, neach);
b1 = linspace(Cstart*sind(hue1),    0, neach);

L2 = linspace(Lend, Lstart            , neach);
a2 = linspace(   0, Cstart*cosd(hue2), neach);
b2 = linspace(   0, Cstart*sind(hue2), neach);

Lab1 = [L1' a1' b1'];
Lab2 = [L2' a2' b2'];
Lab  = [Lab1;Lab2];
cmap = soft_lab2rgb(Lab, use_uplab);

img = repmat(cmap,[1 1 20]);
img = permute(img,[1 3 2]);
figure;
imagesc(img);
axis xy;
title(sprintf('Best seperation distance (%d, %d) = %.2f',hue1,hue2,sep_dist));

fprintf('Best 2norm seperation solution, LCh: %s -> %s -> %s\n',...
    mat2str([Lstart Cstart hue1]),...
    mat2str([Lend 0 0]),...
    mat2str([Lstart Cstart hue2]));
fprintf('Best 2norm seperation solution, Lab: %s -> %s -> %s\n',...
    mat2str([Lstart Cstart*cosd(hue1) Cstart*sind(hue1)]),...
    mat2str([Lend 0 0]),...
    mat2str([Lstart Cstart*cosd(hue2) Cstart*sind(hue2)]));

%% Best chroma seperation

% m = best_sep_m;

[n_best_C_sep,I] = max(all_best_C_sep(:));
m = all_best_C_sep_m(I);
[ih1,ih2] = ind2sub(size(all_best_C_sep),I);

hue1 = hue1_range(ih1);
hue2 = hue2_range(ih2);

gh1 = g.lch(g.lch(:,3)==hue1,:);
gh2 = g.lch(g.lch(:,3)==hue2,:);

jointL = intersect(gh1(:,1),gh2(:,1));
gh1 = gh1(ismember(gh1(:,1),jointL),:);
gh2 = gh2(ismember(gh1(:,1),jointL),:);

jointC = min(gh1(:,2),gh2(:,2));
[Cmax,I_Cmax] = max(jointC);
L_Cmax = gh2(I_Cmax,1);

switch via_black
    case 0
        valid_I = I_Cmax:length(jointC);
    case 1
        valid_I = 1:I_Cmax;
    otherwise
        error('Bad via colour');
end
initL = jointL(valid_I);

% % Plot the blue and red chroma curves to see how well they match
% figure;hold on;set(gca,'Color',[.467 .467 .467]);box on;
% plot(gh1(:,1),gh1(:,2),'b-');
% plot(gh2(:,1),gh2(:,2),'r-');
% % plot(jointL,jointC,'kx');
% title(sprintf('Best C max (%d, %d) = %.2f',hue1,hue2,Cstart));


    % Given a gradient, compute the chroma values for each L
%     C = m*L+k
    % Since L=100 has C=0, curve is offset by
%     k = C-m*L
    if via_black
        initk = 0;
    else
        initk = -m*100;
    end
    initC = m*initL+initk;
    
    % Find min chroma offset needed so all datapoints are in gamut
    k = initk - max(initC - jointC(valid_I));
    
    % Find the intersection of this moved line with C=0
%     C = m*L+k
%     0 = m*Lend+k
%     Lend = -k/m
    Lend = -k/m;
    
    % Now we know the L start value needed to make this gradient work
    
    % Find the furthest we can draw the line before intersecting chroma
    % boundary
    if via_black
        moveL = jointL(jointL>=Lend);
    else
        moveL = jointL(jointL<=Lend);
    end
    moveC = m*moveL+k;
    if via_black
        I_start = find(moveC>jointC(jointL>=Lend),1,'first')-1;
        if I_start<0; continue; end;
    else
        I_start = find(moveC>jointC(jointL<=Lend),1,'last')+1;
        if I_start>length(moveC); continue; end;
    end
    Lstart = jointL(I_start);
    Cstart = moveC(I_start);


% % Plot derivation of map
% figure;hold on;set(gca,'Color',[.467 .467 .467]);box on;
% plot(jointL,jointC,'k-');
% plot(initL,initC,'y-');
% plot(Lend,0,'yo');
% plot(moveL,moveC,'g-');
% plot(Lstart,Cstart,'go');
% title(sprintf('Chosen seperation distance (%d, %d) = %.2f',hue1,hue2,sep_dist));

% Plot the blue and red chroma curves to see how well they match
% WITH derivation of map
figure;hold on;set(gca,'Color',[.467 .467 .467]);box on;
plot(gh1(:,1),gh1(:,2),'b-');
plot(gh2(:,1),gh2(:,2),'r-');
plot(initL,initC,'y-');
plot(Lend,0,'yo');
plot(moveL,moveC,'g-');
plot(Lstart,Cstart,'go');
% plot(jointL,jointC,'kx');
title(sprintf('Chosen seperation distance (%d, %d) = %.2f',hue1,hue2,sep_dist));


% Plot sample colormap
neach = 32;
L1 = linspace(Lstart            , Lend, neach);
a1 = linspace(Cstart*cosd(hue1),    0, neach);
b1 = linspace(Cstart*sind(hue1),    0, neach);

L2 = linspace(Lend, Lstart            , neach);
a2 = linspace(   0, Cstart*cosd(hue2), neach);
b2 = linspace(   0, Cstart*sind(hue2), neach);

Lab1 = [L1' a1' b1'];
Lab2 = [L2' a2' b2'];
Lab  = [Lab1;Lab2];
cmap = soft_lab2rgb(Lab, use_uplab);

img = repmat(cmap,[1 1 20]);
img = permute(img,[1 3 2]);
figure;
imagesc(img);
axis xy;
title(sprintf('Best C max (%d, %d) = %.2f',hue1,hue2,Cstart));

fprintf('Best chroma seperation solution, LCh: %s -> %s -> %s\n',...
    mat2str([Lstart Cstart hue1]),...
    mat2str([Lend 0 0]),...
    mat2str([Lstart Cstart hue2]));
fprintf('Best chroma seperation solution, Lab: %s -> %s -> %s\n',...
    mat2str([Lstart Cstart*cosd(hue1) Cstart*sind(hue1)]),...
    mat2str([Lend 0 0]),...
    mat2str([Lstart Cstart*cosd(hue2) Cstart*sind(hue2)]));


%% Chosen hues, for its optimal seperation

if isempty(handpicked_hue1) || isempty(handpicked_hue2)
    return;
end

% Handpicked from the matrix of values.
hue1 = handpicked_hue1; 
hue2 = handpicked_hue2;

ih1 = find(hue1_range==hue1,1);
ih2 = find(hue2_range==hue2,1);

n_best_sep = all_best_sep(ih1,ih2);
m = all_best_sep_m(ih1,ih2);

gh1 = g.lch(g.lch(:,3)==hue1,:);
gh2 = g.lch(g.lch(:,3)==hue2,:);

jointL = intersect(gh1(:,1),gh2(:,1));
gh1 = gh1(ismember(gh1(:,1),jointL),:);
gh2 = gh2(ismember(gh1(:,1),jointL),:);

jointC = min(gh1(:,2),gh2(:,2));
[Cmax,I_Cmax] = max(jointC);
L_Cmax = gh2(I_Cmax,1);

switch via_black
    case 0
        valid_I = I_Cmax:length(jointC);
    case 1
        valid_I = 1:I_Cmax;
    otherwise
        error('Bad via colour');
end
initL = jointL(valid_I);

% % Plot the blue and red chroma curves to see how well they match
% figure;hold on;set(gca,'Color',[.467 .467 .467]);box on;
% plot(gh1(:,1),gh1(:,2),'b-');
% plot(gh2(:,1),gh2(:,2),'r-');
% % plot(jointL,jointC,'kx');
% title(sprintf('Chosen seperation distance (%d, %d) = %.2f',hue1,hue2,sep_dist));


    % Given a gradient, compute the chroma values for each L
%     C = m*L+k
    % Since L=100 has C=0, curve is offset by
%     k = C-m*L
    if via_black
        initk = 0;
    else
        initk = -m*100;
    end
    initC = m*initL+initk;
    
    % Find min chroma offset needed so all datapoints are in gamut
    k = initk - max(initC - jointC(valid_I));
    
    % Find the intersection of this moved line with C=0
%     C = m*L+k
%     0 = m*Lend+k
%     Lend = -k/m
    Lend = -k/m;
    
    % Now we know the L start value needed to make this gradient work
    
    % Find the furthest we can draw the line before intersecting chroma
    % boundary
    if via_black
        moveL = jointL(jointL>=Lend);
    else
        moveL = jointL(jointL<=Lend);
    end
    moveC = m*moveL+k;
    if via_black
        I_start = find(moveC>jointC(jointL>=Lend),1,'first')-1;
        if I_start<0; continue; end;
    else
        I_start = find(moveC>jointC(jointL<=Lend),1,'last')+1;
        if I_start>length(moveC); continue; end;
    end
    Lstart = jointL(I_start);
    Cstart = moveC(I_start);

    % Define metric of how good this gradient is
    % Euclidian distance from L start value and C max value
    sep_dist = sqrt((Lstart-Lend).^2+Cstart.^2);
    

% % Plot derivation of map
% figure;hold on;set(gca,'Color',[.467 .467 .467]);box on;
% plot(jointL,jointC,'k-');
% plot(initL,initC,'y-');
% plot(Lend,0,'yo');
% plot(moveL,moveC,'g-');
% plot(Lstart,Cstart,'go');
% title(sprintf('Chosen seperation distance (%d, %d) = %.2f',hue1,hue2,sep_dist));

% Plot the blue and red chroma curves to see how well they match
% WITH derivation of map
figure;hold on;set(gca,'Color',[.467 .467 .467]);box on;
plot(gh1(:,1),gh1(:,2),'b-');
plot(gh2(:,1),gh2(:,2),'r-');
plot(initL,initC,'y-');
plot(Lend,0,'yo');
plot(moveL,moveC,'g-');
plot(Lstart,Cstart,'go');
% plot(jointL,jointC,'kx');
title(sprintf('Chosen seperation distance (%d, %d) = %.2f',hue1,hue2,sep_dist));

% Plot sample colormap
neach = 32;
L1 = linspace(Lstart           , Lend, neach);
a1 = linspace(Cstart*cosd(hue1),    0, neach);
b1 = linspace(Cstart*sind(hue1),    0, neach);

L2 = linspace(Lend, Lstart           , neach);
a2 = linspace(   0, Cstart*cosd(hue2), neach);
b2 = linspace(   0, Cstart*sind(hue2), neach);

Lab1 = [L1' a1' b1'];
Lab2 = [L2' a2' b2'];
Lab  = [Lab1;Lab2];
cmap = soft_lab2rgb(Lab, use_uplab);

img = repmat(cmap,[1 1 20]);
img = permute(img,[1 3 2]);
figure;
imagesc(img);
axis xy;
title(sprintf('Chosen seperation distance (%d, %d) = %.2f',hue1,hue2,sep_dist));

fprintf('Chosen hues solution, LCh: %s -> %s -> %s\n',...
    mat2str([Lstart Cstart hue1]),...
    mat2str([Lend 0 0]),...
    mat2str([Lstart Cstart hue2]));
fprintf('Chosen hues solution, Lab: %s -> %s -> %s\n',...
    mat2str([Lstart Cstart*cosd(hue1) Cstart*sind(hue1)]),...
    mat2str([Lend 0 0]),...
    mat2str([Lstart Cstart*cosd(hue2) Cstart*sind(hue2)]));
