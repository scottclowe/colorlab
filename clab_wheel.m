%LABWHEEL A set of colors from a circle in CIELAB color space
%   RGB = LABWHEEL(N,L,C,H0) finds a set of N colors equally spaced in
%   CIELAB space in a circle using CIELCHab co-ordinates. The circle is
%   taken anti-clockwise with constant lightness, L and contant chroma, C,
%   and an initial Hue of H0. Returns RGB, an n-by-3 vector of these colors
%   represented in sRGB space.
%   
%   Scott Lowe, April 2013

function [rgb,params] = clab_wheel(n, L, c, h0, h1, dbg)

% -------------------------------------------------------------------------
if nargin<6
    dbg = true;
end
if nargin<4
    h0 = 0;
end
if nargin<5
    h1 = h0+2*pi;
end
if nargin<3
    c = 60; % chroma
end
if nargin<2
    L = 60; % Lightness
    L = 73.875;
end
if nargin>1 && nargin<4
    % Find from srgb gamut the largest chroma which is available for all H
    rgbgamut = fetch_cielchab_gamut('srgb');
    c = min(rgbgamut.lchmesh.cgrid(:,rgbgamut.lchmesh.Lvec==L));
end

% Default with same number of colors as in use for current colormap
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1);
end

% -------------------------------------------------------------------------
use_uplab = false;

% -------------------------------------------------------------------------
h = linspace(h0, h1, (n+1)).';
h = h(1:end-1);

a = c*cosd(h);
b = c*sind(h);
LL = repmat(L,n,1);

Lab = [LL a b];

rgb = hard_lab2rgb(Lab, use_uplab);

% -------------------------------------------------------------------------
% If dbg mode, display a figure of the outputted colormap
if dbg;
    img = repmat(rgb,[1 1 20]);
    img = permute(img,[1 3 2]);
    figure;
    imagesc(img);
    axis xy;
    title('Output colormap');
    plot_labcurve_rgbgamut(Lab, use_uplab);
end

% -------------------------------------------------------------------------
if nargout<2; return; end
params.n  = n;
params.L  = L;
params.c  = c;
params.h0 = h0;
params.h1 = h1;

end