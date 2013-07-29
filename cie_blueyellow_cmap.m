function cmap = cie_blueyellow_cmap(n, attr, func, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<4 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end
if nargin<3
    func = []; % function to map from cielab to srgb
end
if nargin<2 || isempty(attr)
    attr = []; % Colormap type option
end
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1); % Number of colours in the colormap
end

% -------------------------------------------------------------------------
% CIE  [  L*   a*   b*]
lab1 = [  58   0   -65];
lab2 = [  89   0    87];

% -------------------------------------------------------------------------
L = linspace(lab1(1), lab2(1), n);
a = linspace(lab1(2), lab2(2), n);
b = linspace(lab1(3), lab2(3), n);

Lab = [L' a' b'];

cmap = gd_lab2rgb(Lab, func);

% -------------------------------------------------------------------------
% If dbg mode, display a figure of the outputted colormap
if dbg;
    img = repmat(cmap,[1 1 20]);
    img = permute(img,[1 3 2]);
    figure;
    imagesc(img);
    axis xy;
    title('Output colormap');
end

end