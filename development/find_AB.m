clear all;

%% Parameters

use_uplab = false;

% hue1_range = 260:315; %0:179;
% hue2_range =  10:70; %180+[-15:15];

% hue1_range = 292:.25:305; %260:315; %0:359; %260:315;
% hue2_range =  39:.25:41;  % 10:70;  %0:359; % 10:70 ;

% hue1_range = 307:.25:312; %200:320; %260:315; %0:359; %260:315;
% hue2_range =  52:.25:56; % 20:90;  % 10:70;  %0:359; % 10:70 ;
% h2_isdynamic = 0; % If true then h2=h1+h2

% % Find best separation of two opposing hues
% hue1_range = 0:179;
% hue2_range = 180;
% h2_isdynamic = 1; % If true then h2=h1+h2
% L_range = [];
% % Outcome: purple-green

hue1_range = 80:100;
hue2_range = 260:280;
h2_isdynamic = 0;
L_range = [];


hue1_range =  60:103;
hue2_range = 180:308;
h2_isdynamic = 0;
L_range = [];


hue1_range = [102:.25:103];
hue2_range = 270:307;
h2_isdynamic = 0;
L_range = [];

%% Main code

disp('------------------------------------------------------------------');
disp('CIELab rgb-gamut point seperation maximiser');
disp(['Hue 1 range: ' num2str(hue1_range(1)) ' - ' num2str(hue1_range(end))]);
if ~isempty(hue2_range); disp(['Hue 2 range: ' num2str(hue2_range(1)) ' - ' num2str(hue2_range(end))]); end;
if h2_isdynamic; disp('Hue 2 is defined relative to hue 1'); end;
disp('------------------------------------------------------------------');

hue1_range = mod(hue1_range,360);

g = fetch_cielchab_gamut('srgb', [], [], use_uplab);

all_sep      = nan(length(hue1_range),min(length(hue2_range),1));
all_Cjoint   = nan(length(hue1_range),min(length(hue2_range),1));
all_Cjoint_L = nan(length(hue1_range),min(length(hue2_range),1));
all_Csum     = nan(length(hue1_range),min(length(hue2_range),1));
all_Csum_L   = nan(length(hue1_range),min(length(hue2_range),1));
all_h1       = nan(length(hue1_range),min(length(hue2_range),1));
all_h2       = nan(length(hue1_range),min(length(hue2_range),1));

% Lab euclidian distance
for ih1=1:length(hue1_range)
    hue1 = hue1_range(ih1);
    % If no second hue, we force the second hue to be opposing the first
    if isempty(hue2_range)
        my_hue2_range = mod(hue1+180,360);
    elseif h2_isdynamic
        my_hue2_range = mod(hue1+hue2_range,360);
    else
        my_hue2_range = hue2_range;
    end
    
    for ih2=1:length(my_hue2_range)
        hue2 = my_hue2_range(ih2);
        
        all_h1(ih1,ih2) = hue1;
        all_h2(ih1,ih2) = hue2;
        
        gh1 = g.lch(g.lch(:,3)==hue1,:);
        gh2 = g.lch(g.lch(:,3)==hue2,:);
        
        jointL = intersect(gh1(:,1),gh2(:,1));
        
        if ~isempty(L_range)
%             jointL = intersect(jointL,L_range);
            jointL = jointL(jointL>=min(L_range)&jointL<=max(L_range));
        end
        
        gh1 = gh1(ismember(gh1(:,1),jointL),:);
        gh2 = gh2(ismember(gh2(:,1),jointL),:);
        
        
        
        % 2norm
        
        [C1max,I_C1max] = max(gh1(:,2));
        [C2max,I_C2max] = max(gh2(:,2));
        
        Lch1 = gh1(I_C1max,:);
        Lch2 = gh2(I_C2max,:);
        
        Lab1 = [Lch1(1) Lch1(2)*cosd(Lch1(3)) Lch1(2)*sind(Lch1(3))];
        Lab2 = [Lch2(1) Lch2(2)*cosd(Lch2(3)) Lch2(2)*sind(Lch2(3))];
        
        sep = sqrt(sum((Lab1-Lab2).^2,2));
        
        all_sep(ih1,ih2) = sep;
        
        
        % Sum C
        [C_sum_max,I_Csum_max] = max(gh1(:,2)+gh2(:,2));
        L_Csum = gh2(I_Csum_max,1);
        
        all_Csum(ih1,ih2) = C_sum_max;
        all_Csum_L(ih1,ih2) = L_Csum;
        
        
        % Joint C
        jointC = min(gh1(:,2),gh2(:,2));
        [Cmax,I_Cmax] = max(jointC);
        L_Cmax = gh2(I_Cmax,1);
        
        all_Cjoint(ih1,ih2) = Cmax;
        all_Cjoint_L(ih1,ih2) = L_Cmax;
    end
end

%%

[best_sep,I] = max(all_sep(:));
best_sep_h1 = all_h1(I);
best_sep_h2 = all_h2(I);

[best_Csum,I] = max(all_Csum(:));
best_Csum_h1 = all_h1(I);
best_Csum_h2 = all_h2(I);
best_Csum_L  = all_Csum_L(I);

[best_Cjoint,I] = max(all_Cjoint(:));
best_Cjoint_h1 = all_h1(I);
best_Cjoint_h2 = all_h2(I);
best_Cjoint_L  = all_Cjoint_L(I);

fprintf('Best Euclidian:    sep(%2$3.2f, %3$3.2f)=%1$3.2f\n'           , best_sep,   best_sep_h1,   best_sep_h2 );
fprintf('Best chroma sum:   sep(%2$3.2f, %3$3.2f)=%1$3.2f at L=%4$3.2f\n', best_Csum,  best_Csum_h1,  best_Csum_h2,  best_Csum_L);
fprintf('Best chroma joint: sep(%2$3.2f, %3$3.2f)=%1$3.2f at L=%4$3.2f\n', best_Cjoint,best_Cjoint_h1,best_Cjoint_h2,best_Cjoint_L);

%%

if length(hue1_range)>1 && length(hue2_range)>1
    
    figure;

    subplot(1,3,1);
    imagesc(hue2_range, hue1_range, all_sep);
    axis xy;
    title([{'Best Euclidian'} sprintf('sep(%2$3.2f, %3$3.2f)=%1$3.2f\n', best_sep, best_sep_h1, best_sep_h2 )]);
    colormap(cie_hot);
    colorbar;
    ylabel('Hue 1');
    if h2_isdynamic
        xlabel('Hue 2 (offset from h1)');
    else
        xlabel('Hue 2');
    end

    subplot(1,3,2);
    imagesc(hue2_range, hue1_range, all_Csum);
    axis xy;
    title([{'Best chroma sum'} sprintf('sep(%2$3.2f, %3$3.2f)=%1$3.2f\n', best_Csum,  best_Csum_h1,  best_Csum_h2)]);
    colormap(cie_hot);
    colorbar;
    ylabel('Hue 1');
    if h2_isdynamic
        xlabel('Hue 2 (offset from h1)');
    else
        xlabel('Hue 2');
    end

    subplot(1,3,3);
    imagesc(hue2_range, hue1_range, all_Cjoint);
    axis xy;
    title([{'Best chroma joint'} sprintf('sep(%2$3.2f, %3$3.2f)=%1$3.2f\n', best_Cjoint,best_Cjoint_h1,best_Cjoint_h2)]);
    colormap(cie_hot);
    colorbar;
    ylabel('Hue 1');
    if h2_isdynamic
        xlabel('Hue 2 (offset from h1)');
    else
        xlabel('Hue 2');
    end


else

    if length(hue2_range)<=1
        x = hue1_range;
        xlab = 'Hue 1';
    else
        x = hue2_range;
        xlab = 'Hue 2';
    end

    figure;

    subplot(1,3,1);
    plot(x, all_sep);
    ylabel(sprintf('Best Euclidian: sep(%2$3.2f, %3$3.2f)=%1$3.2f\n', best_sep, best_sep_h1, best_sep_h2 ));
    xlabel(xlab);
    xlim(x([1 end]));

    subplot(1,3,2);
    plot(x, all_Csum);
    ylabel(sprintf('Best chroma sum:   sep(%2$3.2f, %3$3.2f)=%1$3.2f\n', best_Csum,  best_Csum_h1,  best_Csum_h2));
    xlabel(xlab);
    xlim(x([1 end]));

    subplot(1,3,3);
    plot(x, all_Cjoint);
    ylabel(sprintf('Best chroma joint: sep(%2$3.2f, %3$3.2f)=%1$3.2f\n', best_Cjoint,best_Cjoint_h1,best_Cjoint_h2));
    xlabel(xlab);
    xlim(x([1 end]));


end

%%
return;
%% Example of manual extraction of particular hues

h1 = 102;
h2 = 282;

if h2_isdynamic

    ih1 = find(hue1_range==h1);
    ih2 = find(hue2_range==(h2-h1));
    
    C = all_Cjoint(ih1,ih2);
    L = all_Cjoint_L(ih1,ih2);
    
    Lch_start = [L C h1];
    Lch_end   = [L C h2];
    
    fprintf('Chosen h1 start LCh: %s\n',mat2str(Lch_start));
    fprintf('Chosen h2   end LCh: %s\n',mat2str(Lch_end));
    
else

    ih1 = find(hue1_range==h1);
    ih2 = find(hue2_range==h2);
    
    C = all_Cjoint(ih1,ih2);
    L = all_Cjoint_L(ih1,ih2);
    
    Lch_start = [L C h1];
    Lch_end   = [L C h2];
    
    fprintf('Chosen h1 start LCh: %s\n',mat2str(Lch_start));
    fprintf('Chosen h2   end LCh: %s\n',mat2str(Lch_end));
    
end
