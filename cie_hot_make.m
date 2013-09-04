function cmap = cie_hot_make(n, attr, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<3 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end
if nargin<2
    attr = []; % Unused
end
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1); % Number of colours in the colormap
end

% -------------------------------------------------------------------------
% Get params
params = get_params();
params.n = n;

% -------------------------------------------------------------------------
cmap = cmap_pinchedspiral_make(params, dbg);

end

%==========================================================================
function params = get_params()

params.use_uplab = false;
params.typ       = 'pow'; % 'pow' or 'sin'
params.expnt     = 2;     % exponent to use for chroma curve
params.L_off     = 0;     % offset to exponent (not for sine curve)
params.c0        = 0.03;  % proportion of maxc to add to all chroma (maxc\approx75, 67<maxc<77)
params.Lmin      = 4;     % minimum lightness (start)
params.Lmax      = 96;    % maximum lightness (end)
params.h1        =  9;
params.h2        = 95;
params.maxc      = 76.2;

end

%==========================================================================
function params = get_params_old(attr)

% These are based on old gamut and may try to go out of bounds

% Initial values
params.use_uplab = false;
params.typ       = 'pow';
params.expnt     = 2;
params.h1        = 0;
params.L_off     = 0;
params.maxc      = [];
params.Lmin      = 0;
params.Lmax      = 100;
params.c0        = 0;

switch attr
    case '1' %25
        % Based on Set#5
        params.expnt  = 2;
        params.h1     = 6.5;
        params.h2     = 106;
        params.L_off  = 0;
        params.Lmin   = 5; %2;
        params.Lmax   = 95; %98;
        params.maxc   = 73.9; % 75.6 -> 74.25 -> 73.9
        
    case '2' %9
        params.expnt  = 2.2;
        params.h1     = 15;
        params.h2     = 100;
        params.L_off  = 0;
        params.Lmin   = 5; %3;
        params.Lmax   = 95; %97;
        params.maxc   = 73; %[];
        
    case '3' %11
        params.expnt  = 2.3;
        params.h1     = 10;
        params.h2     = 100;
        params.L_off  = 5;
        params.Lmin   = 5; %3;
        params.Lmax   = 95; %97;
        params.maxc   = 71; %[];
        
    case '4' %10
        params.expnt  = 2.5;
        params.h1     = 15;
        params.h2     = 100;
        params.L_off  = 5;
        params.Lmin   = 5; %3;
        params.Lmax   = 95; %97;
        params.maxc   = 67.8; %[];
        
    case '5' %12
        params.expnt  = 2.3;
        params.h1     = 15;
        params.h2     = 100;
        params.L_off  = 5;
        params.Lmin   = 5; %3;
        params.Lmax   = 95; %97;
        params.maxc   = 67.7; %[];
        
    case 'cool:1' %104
        params.expnt  = 2;
        params.h1     = 340;
        params.h2     = 172;
        params.L_off  = -15;
        params.Lmin   = 3;
        params.Lmax   = 85;
        params.maxc   = 70; %[];
        
    otherwise
        error('Unfamiliar attribute');

end

end