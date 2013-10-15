%LABWHEEL A set of colors from a circle in CIELAB color space
%   RGB = LABWHEEL(N,LL,C,H0) finds a set of N colors equally spaced in
%   CIELAB space in a circle using CIELCHab co-ordinates. The circle is
%   taken anti-clockwise with constant lightness, LL and contant chroma, C,
%   and an initial Hue of H0. Returns RGB, an n-by-3 vector of these colors
%   represented in sRGB space.
%   
%   Scott Lowe, April 2013

function [varargout] = clab_wheel(n, LL, h0, h1, dbg)

% -------------------------------------------------------------------------
use_uplab = false;

% -------------------------------------------------------------------------
% Default with same number of colors as in use for current colormap
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1);
end
if nargin<2 || isempty(LL)
    LL = [];
end
if nargin<3 || isempty(h0)
    h0 = 15; % Default start hue with CIELab
    if use_uplab
        h0 = h0+15; % Hues in UPLab are rotated by about 15 degrees, compared with CIELab
    end
end
if nargin<4 || isempty(h1)
    h1 = h0+360; % Default with a full circle
end
if nargin<5
    dbg = false;
end

% -------------------------------------------------------------------------
% Input parsing
clrstr = false;
if ischar(h0); h0=colorstr2h(h0); clrstr=true; end
if ischar(h1); h1=colorstr2h(h1); clrstr=true; end

if clrstr && (h1<h0 || h0==h1)
    h1 = h1+360;
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

if isempty(LL)
    % If no L is given, find the optimal L giving maximum chroma
    cc = min(rgbgamut.lchmesh.cgrid(hli,:),[],1);
    [cc,Lli] = max(cc);
    LL = rgbgamut.lchmesh.Lvec(Lli);
else
    % If L is given, find optimal chroma for each L given
    Lli = ismember(rgbgamut.lchmesh.Lvec,LL);
    cc = min(rgbgamut.lchmesh.cgrid(hli,Lli));
end

% -------------------------------------------------------------------------

if mod(h1-h0,360)>359 || mod(h1-h0,360)<1
    h = linspace(h0, h1, (n+1)).';
    h = h(1:end-1);
else
    h = linspace(h0, h1, n).';
end

aa = bsxfun(@times, cc, cosd(h));
bb = bsxfun(@times, cc, sind(h));

% -------------------------------------------------------------------------
% Now pick out these colours

LL_Lab = cell(length(LL),1);
LL_rgb = cell(length(LL),1);

for j=1:length(LL)
    L = repmat(LL(j),n,1);
    a = aa(:,j);
    b = bb(:,j);
    Lab = [L a b];
    rgb = hard_lab2rgb(Lab, use_uplab);
    LL_Lab{j} = Lab;
    LL_rgb{j} = rgb;
end

% -------------------------------------------------------------------------
% Output all in one matrix, or in individual matrices for each L
% depending on number of outputs
if nargout>1
    varargout = LL_rgb;
else
    rgb = cell2mat(LL_rgb);
    I = bsxfun(@plus,(1:n),n*(0:(length(LL)-1))');
    varargout = {rgb(I(:),:)};
end

% -------------------------------------------------------------------------
% If dbg mode, display a figure of the outputted colormap
if dbg
    rgb = cell2mat(varargout);
    figure;
    imagesc(permute(rgb,[1 3 2]));
    axis xy;
    title('Output colormap');
    
    LL_Lab2 = nan(n*length(LL),3);
    I = setdiff(1:((n+1)*length(LL)),(n+1)*(1:length(LL)));
    LL_Lab2(I,:) = cell2mat(LL_Lab);
    plot_labcurve_rgbgamut(LL_Lab2, use_uplab);
end
% -------------------------------------------------------------------------
% if nargout<2; return; end
% params.n  = n;
% params.LL  = LL;
% params.c  = c;
% params.h0 = h0;
% params.h1 = h1;
% params.use_uplab = use_uplab;

end

% clab_wheel(32, [], 15, 'y', true);