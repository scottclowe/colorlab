function cmap = clab_yellowblue(n, attr, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<3 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end
if nargin<2 || isempty(attr)
    attr = '';
end
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1); % Number of colours in the colormap
end

% -------------------------------------------------------------------------
makefun = [mfilename '_make'];
makefun = str2func(makefun);

cmap = load_cmap(n, makefun, attr, dbg);

end