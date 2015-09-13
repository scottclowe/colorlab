%% Figure comparing with colorbrewer

if ~exist('colorbrewer.mat','file');
    disp('Cannot locate colorbrewer.mat file');
    disp('This is needed to compare clab with cbrewer');
    return;
end

CB = load('colorbrewer.mat');

brewlabpairings = {...
    {'seq' ,'Greys'}    , {'blackwhite'     ,''} ;...
    {'seq' ,'Purples'}  , {'purple'         ,''} ;...
    {'seq' ,'YlOrRd'}   , {'hot'            ,''} ;...
    {'div' ,'RdBu'}     , {'bluewhitered'   ,''} ;...
    {'div' ,'RdYlBu'}   , {'hotcold'        ,''} ;...
    {'div' ,'RdBu'}     , {'hotcold'        ,''} ;...
    {'qual','Pastel2'}  , {'Lfix_disting'   ,''} ;...
    {'qual','Set2'}     , {'Lfix_disting'   ,''} ;...
    {'qual','Dark2'}    , {'Lfix_disting'   ,''} ;...
    {'qual','Set1'}     , {'Lfix_disting'   ,''} ;...
    {'qual','Paired'}   , {'Lfix_disting'   ,''} ...
%     {'qual','Pastel1'}  , {'Lfix_disting'   ,''} ;...
    };


for i=1:size(brewlabpairings,1)
    hf = figure;
    set(hf,'Color',[.4663 .4663 .4663]);
    
    cbtyp = brewlabpairings{i,1}{1};
    cbnam = brewlabpairings{i,1}{2};
    cbrgb = CB.colorbrewer.(cbtyp).(cbnam){end}/255;
    cbrgb = flipud(cbrgb);
    ncols = size(cbrgb,1);
    
    clnam = brewlabpairings{i,2}{1};
    clatr = brewlabpairings{i,2}{2};
    if strcmpi(clnam,'Lfix_disting') && isempty(clatr);
        cblab = rgb2lab(cbrgb);
        if strcmpi(cbnam,'Paired')
            clatr = [mean(cblab(1:2:end,1)) mean(cblab(2:2:end,1))];
            ncols = ncols/2;
        else
            clatr = mean(cblab(:,1));
        end
    end
    clfun = str2func(['clab_' clnam]);
    clrgb = clfun(ncols,clatr);
    
    ax1 = subplot(1,2,1);
    imagesc(permute(cbrgb,[1 3 2]));
    axis off;
    title(['cbrewer: ' cbnam],'Interpreter','none');
    
    ax2 = subplot(1,2,2);
    imagesc(permute(clrgb,[1 3 2]));
    axis off;
    title(['clab: ' clnam ' ' mat2str(clatr)],'Interpreter','none');
    
    set([ax1 ax2],'YDir','normal');
end
