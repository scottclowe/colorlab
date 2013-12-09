function cmap = clab_hotcold_make(n, attr, dbg)

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
% Parameters

% Increment n to be an odd number if it is even
n = n + ~mod(n,2);

% Make a parameter structure to pass to builder
params = struct('n',n);


% -------------------------------------------------------------------------
% Parameter Set #1
% params.h1edg = 316;
% params.h1mid = 270;
% params.h2edg = 37;
% params.h2mid = 80;
% params.Ledg  = 22;
% params.Lmid  = 94;
% params.Lmaxc = 43.5208;
% params.c0    = 0;
% params.maxc  = 70.2095;
% params.typ   = 'sin';
% params.expnt = 1;
% params.use_uplab = false;

% Parameter Set #2: More diverse colours
% params.h1edg = 318;
% params.h1mid = 270;
% params.h2edg = 27;
% params.h2mid = 92;
% params.Ledg  = 17;
% params.Lmid  = 92;
% params.Lmaxc = 44.8764;
% params.c0    = 0;
% params.maxc  = 71.2675;
% params.typ   = 'sin';
% params.expnt = 1;
% params.use_uplab = false;

% Parameter Set #3: Best yet
params.h1edg = 318;
params.h1mid = 270;
params.h2edg = 23;
params.h2mid = 90;
params.Ledg  = 19;
params.Lmid  = 88;
params.Lmaxc = 45.5438;
params.maxc  = 75; %75.0389;
params.c0    = 0;
params.typ   = 'sin';
params.expnt = 1;
params.use_uplab = false;

% -------------------------------------------------------------------------
cmap = makecmap_AwpBtwist(params, dbg);

end