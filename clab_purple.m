function cmap = clab_purple(n, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<2 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1); % Number of colours in the colormap
end

% -------------------------------------------------------------------------
makefun = [mfilename '_make'];
makefun = str2func(makefun);

cmap = load_cmap(n, makefun, [], dbg);

end