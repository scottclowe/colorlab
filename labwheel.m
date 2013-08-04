%LABWHEEL A set of colors from a circle in CIELAB color space
%   RGB = LABWHEEL(N,L,C,H0) finds a set of N colors equally spaced in
%   CIELAB space in a circle using CIELCHab co-ordinates. The circle is
%   taken anti-clockwise with constant lightness, L and contant chroma, C,
%   and an initial Hue of H0. Returns RGB, an n-by-3 vector of these colors
%   represented in sRGB space.
%   
%   Scott Lowe, April 2013

function [rgb,params] = labwheel(n,L,c,h0,h1)

% -------------------------------------------------------------------------
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
end
if nargin>1 && nargin<4
    % Find from srgb gamut the largest chroma which is available for all H
    gamut = fetch_cielchab_gamut('srgb');
    g = gamut.lch(gamut.lch(:,1)==L,:);
    c = min(g(:,2));
end
% if nargin<=1
%     L = 50;
%     c = 55;
%     h0 = 4.6600-2*pi; % 267 deg = 4.6600 rad
%     h1 = 2.5656;      % 147 deg = 2.5656 rad
% end
if nargin<=1
    % If no input, use hand picked set with partial circle of hues
    L = 60;
    c = 60;
    h0 = 4.60; % 264/360*2*pi-2*pi; % 264 deg = 4.6077 rad
    h1 = 2.62; % 150/360*2*pi;      % 150 deg = 2.6180 rad
end
% Default with same number of colors as in use for current colormap
if nargin<1
    n = size(get(gcf,'colormap'),1);
end

% -------------------------------------------------------------------------
use_uplab = false;

% -------------------------------------------------------------------------
h = linspace(h0, h1, (n+1)).';
h = h(1:end-1);

a = c*cos(h);
b = c*sin(h);
Ls = repmat(L,n,1);

Lab = [Ls, a, b];

rgb = gd_lab2rgb(Lab, use_uplab);

% -------------------------------------------------------------------------
if nargout<2; return; end
params.n  = n;
params.L  = L;
params.c  = c;
params.h0 = h0;
params.h1 = h1;

end