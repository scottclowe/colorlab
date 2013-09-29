function cmap = clab_purple_make(n, attr, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<3 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end
if nargin<2
    attr = ''; % Unused
end
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1); % Number of colours in the colormap
end

% -------------------------------------------------------------------------
% Get params
params = get_params();
params.n = n;

% -------------------------------------------------------------------------
cmap = makecmap_pinchedspiral(params, dbg);

end

%==========================================================================
function params = get_params()

params.use_uplab = false;
params.typ       = 'pow'; % 'pow' or 'sin'
params.expnt     = 2;     % exponent to use for chroma curve
params.L_off     = 0;     % offset to exponent (not for sine curve)
params.c0        = 0.025; % proportion of maxc to add to all chroma (maxc\approx75, 67<maxc<77)
params.Lmin      = 4;     % minimum lightness (start)
params.Lmax      = 92;    % maximum lightness (end)
params.h1        = 308.5;
params.h2        = 331;
params.maxc      = 105.9;

end