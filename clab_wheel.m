%LABWHEEL A set of colors from a circle in CIELAB color space
%   RGB = LABWHEEL(N,L,C,H0) finds a set of N colors equally spaced in
%   CIELAB space in a circle using CIELCHab co-ordinates. The circle is
%   taken anti-clockwise with constant lightness, L and contant chroma, C,
%   and an initial Hue of H0. Returns RGB, an n-by-3 vector of these colors
%   represented in sRGB space.
%   
%   Scott Lowe, April 2013

function [rgb,params] = clab_wheel(n, L, h0, h1, dbg)

% -------------------------------------------------------------------------
use_uplab = false;

% -------------------------------------------------------------------------
% Default with same number of colors as in use for current colormap
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1);
end
if nargin<2 || isempty(L)
%     L = 50; % Lightness
%     L = 60; % Lightness
%     L = 73.875;
    L = [];
end
if nargin<3 || isempty(h0)
    h0 = 15;
    if use_uplab
        h0 = h0+15;
    end
end
if nargin<4 || isempty(h1)
    h1 = h0+360;
end
if nargin<5
    dbg = true;
end

% -------------------------------------------------------------------------
rgbgamut = fetch_cielchab_gamut('srgb');

% Find out which h values are included
if h1<h0
    h_srt = h1;
    h_end = h0;
else
    h_srt = h0;
    h_end = h1;
end
if h_end-h_srt >=360
    h_srt = 0;
    h_end = 360;
elseif h_end > 360
    h_srt = h_srt - 360;
    h_end = h_end - 360;
end

% Find the largest chroma which is available for all H
if h_srt>=0
    hli = rgbgamut.lchmesh.hvec >= h_srt & rgbgamut.lchmesh.hvec <= h_end;
else
    hli = rgbgamut.lchmesh.hvec <= h_end | rgbgamut.lchmesh.hvec >= (h_srt+360);
end

if isempty(L)
    cc = min(rgbgamut.lchmesh.cgrid(hli,:),[],1);
    [c,Lli] = max(cc);
    L = rgbgamut.lchmesh.Lvec(Lli);
else
    Lli = rgbgamut.lchmesh.Lvec==L;
    c = min(rgbgamut.lchmesh.cgrid(hli,Lli));
end

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
params.use_uplab = use_uplab;

end