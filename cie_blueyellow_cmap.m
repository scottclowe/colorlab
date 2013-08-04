function cmap = cie_blueyellow_cmap(n, attributes, spacefun, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<4 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end
if nargin<3
    spacefun = []; % function to map from cielab to srgb
end
if nargin<2 || isempty(attributes)
    attributes = 'a0'; % Colormap type option
end
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1); % Number of colours in the colormap
end

% -------------------------------------------------------------------------
switch attributes
    case 'a0'
        % CIE  [  L*    a*  b*]
        lab1 = [  58    0  -65];
        lab2 = [  89    0   87];
    case {'a0:lfix','a0:lfix:1'}
        % NB: The local maxima for joint-chroma is at h=88 (L=59.5, C=63)
        % So this is pretty much optimal for fixed L
        % CIE  [  L*    a*  b*]
        lab1 = [  59.5  0  -62];
        lab2 = [  59.5  0   62];
%     case 'a0:lfix:2'
%         % CIE  [  L*    a*  b*]
%         lab1 = [ 57.75  0  -65];
%         lab2 = [ 57.75  0   61];
%     case 'maxsep'
% %         lch1 = [ 30.5   130    302];
% %         lch2 = [ 90     100    122];
%         lab1 = [ 30.50  68 -109];
%         lab2 = [ 90.00 -52   84];
    case 'primaries'
%         lch1 = [29.5628  131.2104  301.3643];
%         lch2 = [97.6067   94.7101   99.5585];
        lab1 = [ 29.5628   68.2921 -112.0373];
        lab2 = [ 97.6067  -15.7270   93.3952];
    otherwise
        error('Unfamiliar colormap attribute: %s',attributes);
end

% -------------------------------------------------------------------------
L = linspace(lab1(1), lab2(1), n);
a = linspace(lab1(2), lab2(2), n);
b = linspace(lab1(3), lab2(3), n);

Lab = [L' a' b'];

cmap = gd_lab2rgb(Lab, spacefun);

% -------------------------------------------------------------------------
% If dbg mode, display a figure of the outputted colormap
if dbg;
    img = repmat(cmap,[1 1 20]);
    img = permute(img,[1 3 2]);
    figure;
    imagesc(img);
    axis xy;
    title('Output colormap');
    
    plot_labcurve_rgbgamut(Lab)
end

end