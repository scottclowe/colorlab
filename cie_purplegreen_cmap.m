function cmap = cie_purplegreen_cmap(n, attributes, spacefun, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<4 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end
if nargin<3
    spacefun = []; % function to map from cielab to srgb
end
if nargin<2 || isempty(attributes)
    attributes = 'Lfix'; % Colormap type option
end
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1); % Number of colours in the colormap
end

% -------------------------------------------------------------------------
use_uplab = false;

switch lower(attributes)
    case 'lvar'
        % CIE  [  L*      C      h]
        lch1 = [ 46.25  112.5  314]; % purple
        lch2 = [ 87.75  112.5  134]; % green
    case 'lfix'
        % CIE  [  L*     C      h]
        lch1 = [ 60.75   83    314]; % purple
        lch2 = [ 60.75   83    134]; % green
    case 'maxsep'
        lch1 = [ 30.5   130    302]; % blue
        lch2 = [ 90     100    122]; % green
    otherwise
        error('Unfamiliar colormap attribute: %s',attributes);
end

% -------------------------------------------------------------------------
lab1 = [lch1(1) lch1(2).*cosd(lch1(3)) lch1(2).*sind(lch1(3))];
lab2 = [lch2(1) lch2(2).*cosd(lch2(3)) lch2(2).*sind(lch2(3))];

L = linspace(lab1(1), lab2(1), n);
a = linspace(lab1(2), lab2(2), n);
b = linspace(lab1(3), lab2(3), n);

Lab = [L' a' b'];

cmap = gd_lab2rgb(Lab, use_uplab, spacefun);

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

% % hue found using gamut_interactive_app

% % Then Chroma and Lightness found using the following:
% k = min(g.lch(g.lch(:,3)==134,2),g.lch(g.lch(:,3)==134+180,2));
% [C,I] = max(k)
% gg = g.lch(g.lch(:,3)==134,:);
% gg(I,1)
% C