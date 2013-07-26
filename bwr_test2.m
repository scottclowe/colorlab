clear all;
close all;


%% Find seperation distance for a set of different hues

hue_b = 292;
hue_r =  40;

g = fetch_cielchab_gamut('srgb');

hue_b_range = 260:315;
hue_r_range =  10:70 ;

all_best_sep     = nan(length(hue_b_range),length(hue_r_range));
all_best_sep_m   = nan(length(hue_b_range),length(hue_r_range));
all_best_C_sep   = nan(length(hue_b_range),length(hue_r_range));
all_best_C_sep_m = nan(length(hue_b_range),length(hue_r_range));


for ib=1:length(hue_b_range)
    for ir=1:length(hue_r_range)
        hue_b = hue_b_range(ib);
        hue_r = hue_r_range(ir);

ghb = g.lch(g.lch(:,3)==hue_b,:);
ghr = g.lch(g.lch(:,3)==hue_r,:);

jointL = intersect(ghb(:,1),ghr(:,1));
ghb = ghb(ismember(ghb(:,1),jointL),:);
ghr = ghr(ismember(ghb(:,1),jointL),:);

jointC = min(ghb(:,2),ghr(:,2));
[Cmax,I_Cmax] = max(jointC);
L_Cmax = ghr(I_Cmax,1);

% figure;hold on;set(gca,'Color',[.48 .48 .48]);box on;
% plot(ghb(:,1),ghb(:,2),'b-');
% plot(ghr(:,1),ghr(:,2),'r-');
% % plot(jointL,jointC,'kx');


m_segs = diff(jointC(I_Cmax:end))./diff(jointL(I_Cmax:end));
m_max = max(m_segs);
m_min = min(m_segs);

m_list = linspace(m_min,m_max,101);
initL = jointL(I_Cmax:end);

best_sep = 0;
best_sep_m = [];

best_C_sep = 0;
best_C_m = [];

% figure;hold on;set(gca,'Color',[.48 .48 .48]);box on;
% plot(jointL,jointC,'k-');

for i=1:length(m_list)
    m = m_list(i);
    % Given a gradient, compute the chroma values for each L
%     C = m*L+k
    % Since L=100 has C=0, curve is offset by
%     k = C-m*L
    initk = -m*100;
    initC = m*initL+initk;
    
    % Find min chroma offset needed so all datapoints are in gamut
    k = initk - max(initC - jointC(I_Cmax:end));
    
    % Find the intersection of this moved line with C=0
%     C = m*L+k
%     0 = m*Lend+k
%     Lend = -k/m
    Lend = -k/m;
    
    % Now we know the L start value needed to make this gradient work
    
    % Find the furthest we can draw the line before intersecting chroma
    % boundary
    moveL = jointL(jointL<=Lend);
    moveC = m*moveL+k;
    I_start = find(moveC>jointC(jointL<=Lend),1,'last')+1;
    if I_start>length(moveC); continue; end;
    Lstart = jointL(I_start);
    Cstart = moveC(I_start);
    
% %     figure(2);hold on;set(gca,'Color',[.48 .48 .48]);box on;
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


        all_best_sep(ib,ir)     = best_sep;
        all_best_sep_m(ib,ir)   = best_sep_m;
        all_best_C_sep(ib,ir)   = best_C_sep;
        all_best_C_sep_m(ib,ir) = best_C_m;
        
    end
end


%%

figure;
subplot(2,2,1);
imagesc(hue_r_range,hue_b_range,all_best_sep);
colormap(hot); freezeColors; colorbar; cbfreeze;
title('Best seperation distance');
subplot(2,2,2);
imagesc(hue_r_range,hue_b_range,all_best_sep_m);
colormap(hot); freezeColors; colorbar; cbfreeze;
title('gradient');
subplot(2,2,3);
imagesc(hue_r_range,hue_b_range,all_best_C_sep);
colormap(hot); freezeColors; colorbar; cbfreeze;
title('Best chroma distance');
subplot(2,2,4);
imagesc(hue_r_range,hue_b_range,all_best_C_sep_m);
colormap(hot); freezeColors; colorbar; cbfreeze;
title('gradient');

%%

% m = best_sep_m;

[n_best_sep,I] = max(all_best_sep(:));
m = all_best_sep_m(I);
[ib,ir] = ind2sub(size(all_best_sep),I); %numel(all_best_sep));

hue_b = hue_b_range(ib);
hue_r = hue_r_range(ir);

ghb = g.lch(g.lch(:,3)==hue_b,:);
ghr = g.lch(g.lch(:,3)==hue_r,:);

jointL = intersect(ghb(:,1),ghr(:,1));
ghb = ghb(ismember(ghb(:,1),jointL),:);
ghr = ghr(ismember(ghb(:,1),jointL),:);

jointC = min(ghb(:,2),ghr(:,2));
[Cmax,I_Cmax] = max(jointC);
L_Cmax = ghr(I_Cmax,1);
initL = jointL(I_Cmax:end);

% Plot the blue and red chroma curves to see how well they match
figure;hold on;set(gca,'Color',[.48 .48 .48]);box on;
plot(ghb(:,1),ghb(:,2),'b-');
plot(ghr(:,1),ghr(:,2),'r-');
% plot(jointL,jointC,'kx');
title(sprintf('Best seperation distance (%.2f): %d, %d',sep_dist,hue_b,hue_r));

    % Given a gradient, compute the chroma values for each L
%     C = m*L+k
    % Since L=100 has C=0, curve is offset by
%     k = C-m*L
    initk = -m*100;
    initC = m*initL+initk;
    
    % Find min chroma offset needed so all datapoints are in gamut
    k = initk - max(initC - jointC(I_Cmax:end));
    
    % Find the intersection of this moved line with C=0
%     C = m*L+k
%     0 = m*Lend+k
%     Lend = -k/m
    Lend = -k/m;
    
    % Now we know the L start value needed to make this gradient work
    
    % Find the furthest we can draw the line before intersecting chroma
    % boundary
    moveL = jointL(jointL<=Lend);
    moveC = m*moveL+k;
    I_start = find(moveC>jointC(jointL<=Lend),1,'last')+1;
    Lstart = jointL(I_start);
    Cstart = moveC(I_start);

    % Define metric of how good this gradient is
    % Euclidian distance from L start value and C max value
    sep_dist = sqrt((Lstart-Lend).^2+Cstart.^2);
    

% Plot derivation of map
figure;hold on;set(gca,'Color',[.48 .48 .48]);box on;
plot(jointL,jointC,'k-');
plot(initL,initC,'y-');
plot(Lend,0,'yo');
plot(moveL,moveC,'g-');
plot(Lstart,Cstart,'go');
title(sprintf('Best seperation distance (%.2f): %d, %d',sep_dist,hue_b,hue_r));


% Plot sample colormap
neach = 32;
L1 = linspace(Lstart            , Lend, neach);
a1 = linspace(Cstart*cosd(hue_b),    0, neach);
b1 = linspace(Cstart*sind(hue_b),    0, neach);

L2 = linspace(Lend, Lstart            , neach);
a2 = linspace(   0, Cstart*cosd(hue_r), neach);
b2 = linspace(   0, Cstart*sind(hue_r), neach);

Lab1 = [L1' a1' b1'];
Lab2 = [L2' a2' b2'];
Lab  = [Lab1;Lab2];
cform = makecform('lab2srgb');
cmap = applycform(Lab, cform);

img = repmat(cmap,[1 1 20]);
img = permute(img,[1 3 2]);
figure;
imagesc(img);
axis xy;
title(sprintf('Best seperation distance (%.2f): %d, %d',sep_dist,hue_b,hue_r));

%%

% Handpicked from the matrix of values.
hue_b = 289; 
hue_r = 41;

ib = find(hue_b_range==hue_b,1);
ir = find(hue_r_range==hue_r,1);

n_best_sep = all_best_sep(ib,ir);
m = all_best_sep_m(ib,ir);

ghb = g.lch(g.lch(:,3)==hue_b,:);
ghr = g.lch(g.lch(:,3)==hue_r,:);

jointL = intersect(ghb(:,1),ghr(:,1));
ghb = ghb(ismember(ghb(:,1),jointL),:);
ghr = ghr(ismember(ghb(:,1),jointL),:);

jointC = min(ghb(:,2),ghr(:,2));
[Cmax,I_Cmax] = max(jointC);
L_Cmax = ghr(I_Cmax,1);
initL = jointL(I_Cmax:end);

% Plot the blue and red chroma curves to see how well they match
figure;hold on;set(gca,'Color',[.48 .48 .48]);box on;
plot(ghb(:,1),ghb(:,2),'b-');
plot(ghr(:,1),ghr(:,2),'r-');
% plot(jointL,jointC,'kx');
title(sprintf('Chosen seperation distance (%.2f): %d, %d',sep_dist,hue_b,hue_r));


    % Given a gradient, compute the chroma values for each L
%     C = m*L+k
    % Since L=100 has C=0, curve is offset by
%     k = C-m*L
    initk = -m*100;
    initC = m*initL+initk;
    
    % Find min chroma offset needed so all datapoints are in gamut
    k = initk - max(initC - jointC(I_Cmax:end));
    
    % Find the intersection of this moved line with C=0
%     C = m*L+k
%     0 = m*Lend+k
%     Lend = -k/m
    Lend = -k/m;
    
    % Now we know the L start value needed to make this gradient work
    
    % Find the furthest we can draw the line before intersecting chroma
    % boundary
    moveL = jointL(jointL<=Lend);
    moveC = m*moveL+k;
    I_start = find(moveC>jointC(jointL<=Lend),1,'last')+1;
    Lstart = jointL(I_start);
    Cstart = moveC(I_start);

    % Define metric of how good this gradient is
    % Euclidian distance from L start value and C max value
    sep_dist = sqrt((Lstart-Lend).^2+Cstart.^2);
    

% Plot derivation of map
figure;hold on;set(gca,'Color',[.48 .48 .48]);box on;
plot(jointL,jointC,'k-');
plot(initL,initC,'y-');
plot(Lend,0,'yo');
plot(moveL,moveC,'g-');
plot(Lstart,Cstart,'go');
title(sprintf('Chosen seperation distance (%.2f): %d, %d',sep_dist,hue_b,hue_r));


% Plot sample colormap
neach = 32;
L1 = linspace(Lstart            , Lend, neach);
a1 = linspace(Cstart*cosd(hue_b),    0, neach);
b1 = linspace(Cstart*sind(hue_b),    0, neach);

L2 = linspace(Lend, Lstart            , neach);
a2 = linspace(   0, Cstart*cosd(hue_r), neach);
b2 = linspace(   0, Cstart*sind(hue_r), neach);

Lab1 = [L1' a1' b1'];
Lab2 = [L2' a2' b2'];
Lab  = [Lab1;Lab2];
cform = makecform('lab2srgb');
cmap = applycform(Lab, cform);

img = repmat(cmap,[1 1 20]);
img = permute(img,[1 3 2]);
figure;
imagesc(img);
axis xy;
title(sprintf('Chosen seperation distance (%.2f): %d, %d',sep_dist,hue_b,hue_r));

%%

% m = best_sep_m;

[n_best_C_sep,I] = max(all_best_C_sep(:));
m = all_best_C_sep_m(I);
[ib,ir] = ind2sub(size(all_best_C_sep),I);

hue_b = hue_b_range(ib);
hue_r = hue_r_range(ir);

ghb = g.lch(g.lch(:,3)==hue_b,:);
ghr = g.lch(g.lch(:,3)==hue_r,:);

jointL = intersect(ghb(:,1),ghr(:,1));
ghb = ghb(ismember(ghb(:,1),jointL),:);
ghr = ghr(ismember(ghb(:,1),jointL),:);

jointC = min(ghb(:,2),ghr(:,2));
[Cmax,I_Cmax] = max(jointC);
L_Cmax = ghr(I_Cmax,1);
initL = jointL(I_Cmax:end);

% Plot the blue and red chroma curves to see how well they match
figure;hold on;set(gca,'Color',[.48 .48 .48]);box on;
plot(ghb(:,1),ghb(:,2),'b-');
plot(ghr(:,1),ghr(:,2),'r-');
% plot(jointL,jointC,'kx');
title(sprintf('Best C max (%.2f): %d, %d',Cstart,hue_b,hue_r));


    % Given a gradient, compute the chroma values for each L
%     C = m*L+k
    % Since L=100 has C=0, curve is offset by
%     k = C-m*L
    initk = -m*100;
    initC = m*initL+initk;
    
    % Find min chroma offset needed so all datapoints are in gamut
    k = initk - max(initC - jointC(I_Cmax:end));
    
    % Find the intersection of this moved line with C=0
%     C = m*L+k
%     0 = m*Lend+k
%     Lend = -k/m
    Lend = -k/m;
    
    % Now we know the L start value needed to make this gradient work
    
    % Find the furthest we can draw the line before intersecting chroma
    % boundary
    moveL = jointL(jointL<=Lend);
    moveC = m*moveL+k;
    I_start = find(moveC>jointC(jointL<=Lend),1,'last')+1;
    Lstart = jointL(I_start);
    Cstart = moveC(I_start);


% Plot derivation of map
figure;hold on;set(gca,'Color',[.48 .48 .48]);box on;
plot(jointL,jointC,'k-');
plot(initL,initC,'y-');
plot(Lend,0,'yo');
plot(moveL,moveC,'g-');
plot(Lstart,Cstart,'go');
title(sprintf('Best C max (%.2f): %d, %d',Cstart,hue_b,hue_r));


% Plot sample colormap
neach = 32;
L1 = linspace(Lstart            , Lend, neach);
a1 = linspace(Cstart*cosd(hue_b),    0, neach);
b1 = linspace(Cstart*sind(hue_b),    0, neach);

L2 = linspace(Lend, Lstart            , neach);
a2 = linspace(   0, Cstart*cosd(hue_r), neach);
b2 = linspace(   0, Cstart*sind(hue_r), neach);

Lab1 = [L1' a1' b1'];
Lab2 = [L2' a2' b2'];
Lab  = [Lab1;Lab2];
cform = makecform('lab2srgb');
cmap = applycform(Lab, cform);

img = repmat(cmap,[1 1 20]);
img = permute(img,[1 3 2]);
figure;
imagesc(img);
axis xy;
title(sprintf('Best C max (%.2f): %d, %d',Cstart,hue_b,hue_r));
