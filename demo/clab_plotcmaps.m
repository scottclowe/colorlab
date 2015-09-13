
do_save = true;
savefld = fileparts(mfilename('fullpath'));
if ~isempty(savefld); savefld = fullfile(savefld,'figs'); end;

%% Figure showing overview of all colormaps
n0 =  8;
n1 = 32;
n2 = 33;

cmap_settings = {...
    'blackwhite'    , n1, '';
    'purple'        , n1, '';
    'hot'           , n1, '';
    'hotcold'       , n2, '';
    'bluewhitered'  , n2, '';
    'bluewhitered'  , n2, 'clin';
    'bluewhitered'  , n2, 'lfix';
    'greenpurple'   , n2, '';
    'yellowblue'    , n2, '';
    'wheel'         , n0, '';
    'Larc'          , n0, '';
    'Lfix_disting'  , n0, '';
    };

nMap = size(cmap_settings,1);

hf = figure;
set(hf,'Color',[.4663 .4663 .4663]);

iax = 0;
ax1unit = 'normalized'; % [ inches | centimeters | {normalized} | points | pixels | characters ]
ax1xsrt = 0.1;
ax1ysrt = 0.1;
ax1xlen = (1-2*ax1xsrt)/(2*nMap-1);
ax1ylen = 0.8;
ax1xint = ax1xlen;
ax1yint = 0;

ax1 = zeros(nMap, 1);
for iMap=1:nMap
    iax = iax+1;
    ax1(iax) = axes;
    set(ax1(iax),...
        'Units',ax1unit,...
        'Position',[ax1xsrt+(ax1xlen+ax1xint)*(iax-1) ax1ysrt ax1xlen ax1ylen]);
    
    cmapfn = str2func(['clab_' cmap_settings{iMap,1}]);
    cmap = cmapfn(cmap_settings{iMap,2}, cmap_settings{iMap,3});
    
    imagesc(permute(cmap,[1 3 2]));
    set(ax1(iax),...
        'YDir','normal',...
        'XTick',[],...
        'YTick',[]);
    
    axis off;
    
    if isempty(cmap_settings{iMap,3})
        addstr = '';
    else
        addstr = [', ' cmap_settings{iMap,3}];
    end
    hl = ylabel([cmap_settings{iMap,1} ' ( ' num2str(cmap_settings{iMap,2})  addstr ' )'],...
        'Interpreter', 'none',...
        'Visible', 'on');
end

if do_save
    if ~exist(savefld,'dir'); mkdir(savefld); end;
    fnm = fullfile(savefld,'clab_overview');
    if exist('export_fig.m','file')
        export_fig(fnm, '-png', '-eps', '-painters');
    else
        set(gcf, 'InvertHardCopy','off');
        print(fnm,'-dpng','-painters');
    end
end


%% Figure showing LL types

% Linear luminance
% Isoluminant


