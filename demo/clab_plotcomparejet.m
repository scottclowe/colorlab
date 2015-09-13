
do_save = true;
savefld = fileparts(mfilename('fullpath'));
savefld = fullfile(savefld,'figs');

%% Figure comparing jet and clab_hot as a function of number of colours
% Portrait

cmap_settings = {...
    'jet'         ;
    'clab_hot'    ;
    };

nMap = size(cmap_settings,1);

Ncl  = [8 10 12 16 20 24 32 48 64 256];
nNcl = length(Ncl);

hf = figure;
set(hf,'Color',[.4663 .4663 .4663]);

iax = 0;
ax1unit = 'normalized'; % [ inches | centimeters | {normalized} | points | pixels | characters ]
ax1xsrt = 0.1;
ax1ysrt = 0.1;
ax1xint = ax1xsrt;
ax1yint = ax1ysrt;
ax1xlen = (1-2*ax1xsrt+ax1xint)/nNcl;
ax1ylen = (1-2*ax1ysrt+ax1yint)/nMap;
ax1xint = min(ax1xint,ax1xlen/2);
ax1yint = min(ax1yint,ax1ylen/2);
ax1xlen = ax1xlen-ax1xint;
ax1ylen = ax1ylen-ax1yint;

ax1 = zeros(nMap, nNcl);
hyl = zeros(nMap, 1   );
htt = zeros(1   , nNcl);
for iNcl=1:length(Ncl)
    for iMap=1:nMap
        ax1(iMap,iNcl) = axes;
        axpos = [...
            ax1xsrt+(ax1xlen+ax1xint)*(iNcl-1)...
            ax1ysrt+(ax1ylen+ax1yint)*(nMap-iMap)...
            ax1xlen...
            ax1ylen];
        set(ax1(iMap,iNcl),...
            'Units',ax1unit,...
            'Position',axpos);
        
        cmapfn = str2func(cmap_settings{iMap,1});
        cmap = cmapfn(Ncl(iNcl));
        
        imagesc(permute(cmap,[1 3 2]));
        set(ax1(iMap,iNcl),...
            'YDir','normal',...
            'XTick',[],...
            'YTick',[]);
        
        axis off;
        
        if iNcl==1
            hyl(iMap,iNcl) = ylabel(cmap_settings{iMap,1},...
                'Interpreter', 'none',...
                'Visible', 'on');
        end
        if iMap==1
            htt(iMap,iNcl) = title(num2str(Ncl(iNcl)),...
                'Interpreter', 'none',...
                'Visible', 'on');
        end
    end
end

set([htt(:);hyl(:)],...
    'Color',[.1 .1 .1],...
    'Interpreter','none');

if do_save
    if ~exist(savefld,'dir'); mkdir(savefld); end;
    fnm = fullfile(savefld,'compare_ncols_landscape');
    if exist('export_fig.m','file')
        export_fig(fnm, '-png', '-eps', '-painters', '-transparent');
    else
        set(gcf, 'InvertHardCopy','off');
        print(fnm,'-dpng','-painters');
    end
end


%% Figure comparing jet and clab_hot as a function of number of colours
% Landscape

cmap_settings = {...
    'jet'         ;
    'clab_hot'    ;
    };

nMap = size(cmap_settings,1);

Ncl  = [8 10 12 16 20 24 32 48 64 256];
nNcl = length(Ncl);

hf = figure;
set(hf,'Color',[.4663 .4663 .4663]);

iax = 0;
ax1unit = 'normalized'; % [ inches | centimeters | {normalized} | points | pixels | characters ]
ax1xsrt = 0.1;
ax1ysrt = 0.1;
ax1xint = ax1xsrt;
ax1yint = ax1ysrt;
ax1xlen = (1-2*ax1xsrt+ax1xint)/nMap;
ax1ylen = (1-2*ax1ysrt+ax1yint)/nNcl;
ax1xint = min(ax1xint,ax1xlen/2);
ax1yint = min(ax1yint,ax1ylen/2);
ax1xlen = ax1xlen-ax1xint;
ax1ylen = ax1ylen-ax1yint;

ax1 = zeros(nMap, nNcl);
htt = zeros(nMap, 1   );
hyl = zeros(1   , nNcl);
for iNcl=1:nNcl
    for iMap=1:nMap
        ax1(iMap,iNcl) = axes;
        axpos = [...
            ax1xsrt+(ax1xlen+ax1xint)*(iMap-1)...
            ax1ysrt+(ax1ylen+ax1yint)*(nNcl-iNcl)...
            ax1xlen...
            ax1ylen];
        set(ax1(iMap,iNcl),...
            'Units',ax1unit,...
            'Position',axpos);
        
        cmapfn = str2func(cmap_settings{iMap,1});
        cmap = cmapfn(Ncl(iNcl));
        
        imagesc(permute(cmap,[3 1 2]));
        set(ax1(iMap,iNcl),...
            'YDir','normal',...
            'XTick',[],...
            'YTick',[]);
        
        axis off;
        
        if iMap==1
            hyl(iMap,iNcl) = ylabel(num2str(Ncl(iNcl)),...
                'Interpreter', 'none',...
                'Visible', 'on');
        end
        if iNcl==1
            htt(iMap,iNcl) = title(cmap_settings{iMap,1},...
                'Interpreter', 'none',...
                'Visible', 'on');
        end
    end
end

set([htt(:);hyl(:)],...
    'Color',[.1 .1 .1],...
    'Interpreter','none');

if do_save
    if ~exist(savefld,'dir'); mkdir(savefld); end;
    fnm = fullfile(savefld,'compare_ncols_portrait');
    if exist('export_fig.m','file')
        export_fig(fnm, '-png', '-eps', '-painters', '-transparent');
    else
        set(gcf, 'InvertHardCopy','off');
        print(fnm,'-dpng','-painters');
    end
end

%% Figure showing L,C,h,dE for different colormaps

n = 32;

cmap_settings = {...
    'jet'         ;
    'clab_hot'    ;
    };

nMap = size(cmap_settings,1);

% Use the hues for r, g and b primaries, but set to matched L and C
rgbclrs = hard_lab2rgb(lch2lab(clab_cmax(50,[colorstr2h('r') colorstr2h('g') colorstr2h('b')])));

% Show Lightness in grey, a in magenta and b in yellow
labclrs = hard_lab2rgb(lch2lab(...
    vertcat([50 0 0],clab_cmax(50,[0 90]))...
    ));

% Unsure what to use
% Max chroma point   [32.3750 133.5938 306.25]
% Max chroma at L=50 [50      118.5    319   ]
% Min chroma edge    [79.875   43.5    225.5 ]
% Min chroma at L=50 [50       29.4688 214.25]
lchclrs = hard_lab2rgb(lch2lab(...
    vertcat(...
        [33 0 0],...
        [50 118.5 319],...
        [66 0 0])...
    ));

for iMap=1:nMap
    figure;
    
    cmapfn = str2func(cmap_settings{iMap,1});
    cmap = cmapfn(n);
    
    rgb = cmap;
    lab = rgb2lab(rgb);
    lch = lab2lch(lab);
    
    Lab_dif = diff(lab,1,1);
    Lab_sep = sqrt(sum(Lab_dif.^2,2));
    dE1931 = Lab_sep;
    
    dE2000 = ciede2000(lab(2:end,:),lab(1:end-1,:));
    
    subplot(1,4,1);
    hold on;
    for i=1:3
        plot(rgb(:,i),'-','Color',rgbclrs(i,:),'LineWidth',3);
    end
    xlim([0 n+1])
    title('RGB');
    legend('Red','Green','Blue',...
        'Location','SouthOutside');
    
    subplot(1,4,2);
    hold on;
    for i=1:3
        plot(lab(:,i),'-','Color',labclrs(i,:),'LineWidth',3);
    end
    xlim([0 n+1])
    title('L*a*b*');
    legend('Lightness','a*','b*',...
        'Location','SouthOutside');
    
    subplot(1,4,3);
    hold on;
    for i=1:3
        plot(lch(:,i),'-','Color',lchclrs(i,:),'LineWidth',3);
    end
    xlim([0 n+1])
    title('LCh');
    legend('Lightness','Chroma','Hue',...
        'Location','SouthOutside');
    
    subplot(1,4,4);
    hold on;
    plot(dE1931,'-','Color',[.3 .3 .3],'LineWidth',3);
    plot(dE2000,'-','Color',[.7 .7 .7],'LineWidth',3);
    xlim([0 n+1])
    ylim([0 30]);
    title('\DeltaE');
    legend(['Color Dif' 10 '(CIE31)'],['Color Dif' 10 '(CIE2000)'],...
        'Location','SouthOutside');
    
    if do_save
        if ~exist(savefld,'dir'); mkdir(savefld); end;
        fil = sprintf('%s_LCH',cmap_settings{iMap,1});
        fnm = fullfile(savefld,fil);
        if exist('export_fig.m','file')
            export_fig(fnm, '-png', '-eps', '-painters', '-transparent');
        else
            set(gcf, 'InvertHardCopy','off');
            print(fnm,'-dpng','-painters');
        end
    end
end


%% Figures showing curves against the sRGB gamut

n = 32;

cmap_settings = {...
    'jet'         , [], '';
    'clab_hot'    , [], '';
    };

nMap = size(cmap_settings,1);

for iMap=1:nMap
    cmapfn = str2func(cmap_settings{iMap,1});
    if isempty(cmap_settings{iMap,3})
        cmap = cmapfn(n);
    else
        cmap = cmapfn(n, cmap_settings{iMap,3});
    end
    
    rgb = cmap;
    lab = rgb2lab(rgb);
%     lch = lab2lch(lab);
    
    plot_labcurve_rgbgamut(lab, false)
    
    if do_save
        if ~exist(savefld,'dir'); mkdir(savefld); end;
        fil = sprintf('%s_LCH',cmap_settings{iMap,1});
        fnm = fullfile(savefld,fil);
        if exist('export_fig.m','file')
            export_fig(gcf, fnm, '-png', '-eps', '-painters', '-transparent');
        else
            set(gcf, 'InvertHardCopy','off');
            print(fnm,'-dpng','-painters');
        end
    end
    
end

%% Figure showing peaks

dats = cell(1,4);
dat_names = cell(1,4);

load flujet
dats{1} = X;
dat_names{1} = 'flujet';
dats{2} = X - 32;
dat_names{2} = 'flujet-offset';

figure;
image(X);
colormap(jet);

load spine
dats{3} = X;
dat_names{3} = 'spine';

figure;
image(X);
colormap bone;

dats{4}      = peaks(200);
dat_names{4} = 'peaks';

dats{5}      = peaks(200) + randn(200);
dat_names{5} = 'peaks+gaunoise';

dats{6}      = peaks(200) + 2*(rand(200)-0.5);
dat_names{6} = 'peaks+uninoise';


for idat=1:numel(dats)
    
    dat = dats{idat};
    
hf = figure;
set(hf,'Color',[.8 .8 .8]);

Ncl = [8 12 16 32 256];

ax2 = zeros(length(Ncl),4);
% hcb = zeros(length(Ncl),4);
hyl = zeros(length(Ncl),1);
htt = zeros(1          ,4);
nrow = length(Ncl);
ncol = 4;

for iNcl=1:length(Ncl)
    
    ax2(iNcl,1) = subplot(nrow,ncol,1+(iNcl-1)*ncol);
    imagesc(dat);
    colormap(jet(Ncl(iNcl)));
    hcb = colorbar;
    set(hcb,...
        'TickDir','out',...
        'XColor',[.4663 .4663 .4663],...
        'YColor',[.4663 .4663 .4663]);
    freezeColors; cbfreeze;
    if iNcl==1; htt(iNcl,1) = title('MATLAB jet'); end;
    hyl(iNcl,1) = ylabel(['N = ' num2str(Ncl(iNcl))]);

    ax2(iNcl,2) = subplot(nrow,ncol,2+(iNcl-1)*ncol);
    imagesc(dat);
    colormap(hot(Ncl(iNcl)));
    hcb = colorbar;
    set(hcb,...
        'TickDir','out',...
        'XColor',[.4663 .4663 .4663],...
        'YColor',[.4663 .4663 .4663]);
    freezeColors; cbfreeze;
    if iNcl==1; htt(iNcl,2) = title('MATLAB hot'); end;

    ax2(iNcl,3) = subplot(nrow,ncol,3+(iNcl-1)*ncol);
    imagesc(dat);
    colormap(clab_hot(Ncl(iNcl)));
    hcb = colorbar;
    set(hcb,...
        'TickDir','out',...
        'XColor',[.4663 .4663 .4663],...
        'YColor',[.4663 .4663 .4663]);
    freezeColors; cbfreeze;
    if iNcl==1; htt(iNcl,3) = title('clab_hot'); end;

    ax2(iNcl,4) = subplot(nrow,ncol,4+(iNcl-1)*ncol);
    imagesc(dat);
    caxis([-1 1].*max(abs(dat(:))));
    colormap(clab_hotcold(Ncl(iNcl)+~mod(Ncl(iNcl),2))); % Needs to be odd
    hcb = colorbar;
    set(hcb,...
        'TickDir','out',...
        'XColor',[.4663 .4663 .4663],...
        'YColor',[.4663 .4663 .4663]);
    freezeColors; cbfreeze;
    if iNcl==1; htt(iNcl,4) = title('clab_hotcold'); end;

end

set(ax2,...
    'XTick',[],...
    'YTick',[],...
    'XColor',[.4663 .4663 .4663],...
    'YColor',[.4663 .4663 .4663]);
set([htt(:);hyl(:)],...
    'Color',[.1 .1 .1],...
    'Interpreter','none');
linkaxes(ax2);

if do_save
    if ~exist(savefld,'dir'); mkdir(savefld); end;
    fil = sprintf('compare_%s',dat_names{idat});
    fnm = fullfile(savefld,fil);
    if exist('export_fig.m','file')
        export_fig(fnm, '-png', '-eps', '-painters', '-transparent');
    else
        set(gcf, 'InvertHardCopy','off');
        print(fnm,'-dpng','-painters');
    end
end

end
