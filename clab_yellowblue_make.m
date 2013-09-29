function cmap = clab_yellowblue_make(n, attr, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<3 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end
if nargin<2 || isempty(attr)
    attr = ''; % Colormap type option
end
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1); % Number of colours in the colormap
end

% -------------------------------------------------------------------------
% Lookup parameters
params = struct('n',n);
params.use_uplab = false;

switch attr
    case ''
        % CIE  [  L*     C     h ]
%         params.lch1 = [59.375 63.906  98];
%         params.lch2 = [59.375 63.906 278];
        
%         params.lch1 = [64.375 68.8125 102.75];
%         params.lch2 = [64.375 68.8125 307   ];
        
%         params.lch1 = [59.125 64.375 103];
%         params.lch2 = [59.125 64.375 278];
        
        params.lch1 = [59.5   64.5   102];
        params.lch2 = [59.5   64.5   282];

    case {'a0','a0:lfix','a0:lfix:1'}
        % CIE  [  L*     C     h ]
        params.lch1 = [63.375 56.968  90];
        params.lch2 = [63.375 56.968 270];
        
% OLD
%     case 'a0'
%         % CIE  [  L*    a*  b*]
%         lab1 = [  58    0  -65];
%         lab2 = [  89    0   87];
%         
%     case {'a0:lfix','a0:lfix:1'}
%         % NB: The local maxima for joint-chroma is at h=88 (L=59.5, C=63)
%         % So this is pretty much optimal for fixed L
%         % It used to look like a local maxima, but it sure doesn't anymore
%         % CIE  [  L*    a*  b*]
%         lab1 = [  59.5  0  -62];
%         lab2 = [  59.5  0   62];

% OLDER
%     case 'a0:lfix:2'
%         % CIE  [  L*    a*  b*]
%         lab1 = [ 57.75  0  -65];
%         lab2 = [ 57.75  0   61];
%
%     case 'maxsep'
% %         params.lch1 = [ 30.5   130    302];
% %         params.lch2 = [ 90     100    122];
%         lab1 = [ 30.50  68 -109];
%         lab2 = [ 90.00 -52   84];

    case 'primaries'
        
%         lab1 = [97.1279  -21.5484   94.4728]; % yellow
%         lab2 = [32.3006   79.1839 -107.8569]; % blue
        
        params.lch1 = [97.1279   96.8991  102.8489]; % yellow
        params.lch2 = [32.3006  133.8029  306.2845]; % blue
        
    otherwise
        error('Unfamiliar colormap attribute: %s',attr);
end

% -------------------------------------------------------------------------
% Build the colormap
cmap = makecmap_ABlin(params, dbg);

end