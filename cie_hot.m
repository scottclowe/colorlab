function cmap = cie_hot(n, attr, spacefun, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<4 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end
if nargin<3
    spacefun = []; % function to map from cielab to srgb
end
if nargin<2 || isempty(attr)
    attr = '1'; % Colormap type option
end
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1); % Number of colours in the colormap
end

% -------------------------------------------------------------------------
makefun = [mfilename '_make'];
makefun = str2func(makefun);

cmap = load_cmap(n, makefun, attr, spacefun, dbg);

end