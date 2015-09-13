%LABWHEEL A set of colors from a circle in CIELAB color space
%   RGB = CLAB_WHEEL(N,L) finds a set of N colors equally spaced in
%   CIELAB space in a circle using CIELCHab co-ordinates. The circle is
%   taken anti-clockwise with constant lightness of L, and contant chroma
%   which is the maximum possible such that all colours lie within the RGB
%   gamut. The start and end colour is red.
%   red->orange->yellow->green->cyan->blue->purple->magenta->red
%   Output is in sRGB colour space.
%   
%   CLAB_WHEEL with no inputs or N=[] will use the same number of colours
%   as in the current colormap.
%   
%   CLAB_WHEEL(N) or with L=[] will use L as whichever Lightness maximises
%   the chroma of the colours.
%   
%   RGB = CLAB_WHEEL(N,LL) with LL as a vector finds the circle/arc with
%   maximum chroma at each lightness in LL, then interweaves the results to
%   create an N*numel(LL)x3 matrix of colours. The output will go through
%   each lightness with the same hue before progressing to the next hue.
%   
%   [RGB1,RGB2,...] = CLAB_WHEEL(N,LL) with the same number of outputs as
%   there in elements in the vector LL will return the result at each
%   lightness in the corresponding ouptut instead of interspersing them.
%   
%   CLAB_WHEEL(N,LL,H0) with H0 as a scalar sets the start and finish
%   colour to be H0. This can be a colour string such as 'r','y','c'.
%   
%   CLAB_WHEEL(N,LL,H0,H1) will return an arc from the circle which runs
%   from hue H0 to hue H1. Again, H1 may be a colour string. If H1 is more
%   than 360 degrees from H0, more than a complete circle will be returned.
%   If H0 and H1 are both numeric, the direction traversed will respect the
%   values given.
%   If one or both of H0 and H1 are colorstrings, the arc will be taken
%   anti-clockwise from H0 to H1, regardless of the distance between the
%   hues. If you want a clockwise arc, apply flipud to the output.
%   
%   CLAB_WHEEL(N,LL,H) with H as a two-element numeric vector is the same
%   as using CLAB_WHEEL(N,LL,H(1),H(2)).
%   
%   CLAB_WHEEL(...,dEmode) will use dEmode as the colour difference
%   specifications. Available colour distance metrics are:
% 
%      'cie1976' - Uses constant separation in CIELAB space, which is the
%                  same as the CIE1976 dE specifications.
%      'cie2000' - Uses the improved CIE2000 dE specifications, which are
%                  more perceptually uniform.
%      'uplab'   - Uniformly Perceptual LAB space, which maps the Munsell
%                  colours onto CIELAB, allowing a distance metric in
%                  Munsell colours.
%
%   By default, 'cie1976' is used if one more than one Lightness is
%   specified, otherwise 'cie2000' is used.
%   
%   CLAB_WHEEL(...,'matchC') will use the same chroma for all the Lightness
%   values requested with LL. The chroma used will be the minimum chroma
%   from all Lightness used. (Does nothing if LL is scalar or empty.)
%   
%   CLAB_WHEEL(...,'debug') will plot the colours outputted and
%   their locations on the sRGB gamut displayed in the colorspace used
%   (either CIELAB or UPLAB).
%   
%   Examples:
%      cmap = clab_wheel;    % Default settings (Lightly coloured circle)
%      
%      clab_wheel(256,65); % Colours are less soft
%      clab_wheel(256,60,'uplab'); % Colours are less soft
%      
%      % View some more examples. Remove 'debug' in actual use.
%      clab_wheel(10, [55 70], 'debug'); % 10 light and 10 dark shades
%      clab_wheel(10, [55 70], 'cmatch', 'debug'); % same C throughout
%      clab_wheel(256, [], [0 720], 'cie1976', 'debug'); % Does two cycles
%      clab_wheel([], [], 10, 'y', 'debug'); % A "summer" colormap
%      
%      colormap(clab_wheel(...)); % Sets the current colormap

%   Scott Lowe, April 2013

function [varargout] = clab_wheel(n, LL, h0, h1, flg1, flg2, flg3)

% Input Handling ==========================================================
% -------------------------------------------------------------------------

error(nargchk(0, 6, nargin, 'struct'));
vargin2 = {};

% Fixed order inputs ------------------------------------------------------

% Default with same number of colors as in use for current colormap
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1);
end
% Default with chroma-optimal Lightness
if nargin<2 || isempty(LL)
    LL = [];
end

% Fixed input type check
if ~isnumeric(n) && ~isinteger(n)
    error('Number of colours must be an integer');
end
if ischar(LL) % We'll allow clab_wheel(n,'debug');
    vargin2 = {LL};
    LL = [];
end
if ~isnumeric(LL)
    error('Lightness values must be integer');
end
if any(LL<0)
    error('Lightness values must be positive');
end
if any(LL>100)
    error('Lightness values must be between 0 and 100');
end

% Dynamic order inputs ====================================================
% We make a list of inputs so we can prompt the user as to what we expect
% But really, we take the inputs in any order!

% Set varargin variable ---------------------------------------------------
if nargin<3; h0 = []; end;
if nargin<4; h1 = []; end;
if nargin<5; flg1 = []; end;
if nargin<6; flg2 = []; end;
if nargin<7; flg3 = []; end;

vargin2 = {vargin2{:}; h0; h1; flg1; flg2; flg3};

% Set default values ------------------------------------------------------
h0        = [];          % Start hue (default depends on use_uplab)
h1        = [];          % End hue (default depends on h0)
use_uplab = false;       % Whether to use CIELAB not UPLAB (parametrised Munsell)
cie2000   = numel(LL)<2; % Whether to use the CIE2000 colour difference metric
cmatch    = false;       % Whether to use the same Chroma for all Lightness
dbg       = false;       % Whether to debug and plot output colours

% Consider each input in turn ---------------------------------------------
nArgNum = 0; % Number of numeric arguments
for iArg = 1:numel(vargin2)
    arg = vargin2{iArg};
    % Handle numeric inputs and colorstring inputs
    if isnumeric(arg) || ...
            (ischar(arg) && length(arg)==1 && ismember(arg,'rgbwcmyk'))
        nArgNum = nArgNum+1;
        switch nArgNum
            case 1
                h0 = arg;
            case 2
                h1 = arg;
            otherwise
                if isempty(arg); continue; end;
                error('%d is too many numeric and colorstring inputs',nArgNum);
        end
        
    % Handle character string flags
    elseif ischar(arg)
        % Flags
        switch lower(arg)
            case {'dbg','debug','plot','preview','show'}
                dbg = true;
            case {'cmatch','matchc','fixc','cfix'}
                cmatch = true;
            case {'cie1976','cie76'}
                cie2000   = false;
                use_uplab = false;
            case {'cie2000','cie00'}
                cie2000   = true;
                use_uplab = false;
            case {'uplab','use_uplab','munsell'}
                cie2000   = false;
                use_uplab = true;
            otherwise
                error('Unfamiliar flag: %s',arg);
        end
    
    % Otherwise we don't know this class
    else
        error('%s',class(arg));
    end
end

if cie2000 && use_uplab
    error('CIE2000 and UPLAB are both colour difference metrics. Cannot combine them.')
end
if isnumeric(h0) && numel(h0)==2
    h1 = h0(2);
    h0 = h0(1);
end

% Parse dynamic defaults --------------------------------------------------
if isempty(h0)
    h0 = colorstr2h('r',use_uplab); % Start with red
end
if isempty(h1)
    h1 = h0+360; % Default with a full circle
end

% Parse colorstring -------------------------------------------------------
clrstr = false;
if ischar(h0); h0=colorstr2h(h0,use_uplab); clrstr=true; end
if ischar(h1); h1=colorstr2h(h1,use_uplab); clrstr=true; end

% Ensure we go clockwise if hues are colorstrings
if clrstr && (h1<h0 || h0==h1)
    h1 = h1+360;
end

% Main script =============================================================
% -------------------------------------------------------------------------

% Get sRGB gamut boundary
rgbgamut = fetch_cielchab_gamut('srgb',[],[],use_uplab);

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
else
    k = floor(h_end./360);
    h_srt = h_srt - k*360;
    h_end = h_end - k*360;
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
    % If L is set, find optimal chroma for each L given
% %     LL = round(LL/rgbgamut.Lintv)*rgbgamut.Lintv; % Round to Linteveral
    if all(ismember(LL,rgbgamut.lchmesh.Lvec))
        % Method #1 only works with L in the gamut intervals
        Lli = ismember(rgbgamut.lchmesh.Lvec,LL);
        cc = min(rgbgamut.lchmesh.cgrid(hli,Lli));
    else
        % Method #2 
        % Use isingamut to do interpolation
        % It interpolates on the L-h grid to find the maximum C in gamut at each
        % point.
        nBtsHue = 1+ceil((h_end-h_srt)/min(1,rgbgamut.hintv)); % Number of hue bootstraps
        cc = nan(1,numel(LL));
        for iLL=1:numel(LL)
            lch = ones(nBtsHue,1);
            lch(:,1) = LL(iLL);
            lch(:,2) = 1;
            lch(:,3) = linspace(h_srt,h_end,nBtsHue)';
            [TF,lch] = isingamut(lch,rgbgamut,'lch');
            cc(iLL) = min(lch(:,2));
        end
    end
    % Keep the same chroma for all lightness, if needed
    if cmatch
        cc = ones(size(cc))*min(cc(:));
    end
end

% -------------------------------------------------------------------------
% Pick our hues

% If its a whole number of circumnavigations (ish) then we don't want the
% last colour to be the same as the first colour, we want them to be spaced
% out just like the others
if mod(h1-h0,360)>359 || mod(h1-h0,360)<1
    is_circ = 1;
else
    is_circ = 0;
end

if cie2000
    % Bootstrap the arc so we can smooth dE separation
    nBtsDE = ceil(max((h_end-h_srt)*2, n*20));
    h = nan(n+is_circ, numel(LL));
    % Have to do each Lightness separately
    for iLL=1:numel(LL)
        lch = nan(nBtsDE,3);
        lch(:,1) = LL(iLL);
        lch(:,2) = cc(iLL);
        lch(:,3) = linspace(h0, h1, nBtsDE);
        lch = clab_smoothDE(lch, n+is_circ, 'lch', 'cie2000');
        h(:,iLL) = lch(:,3);
    end
    % All the hues should be the same
    h = mean(h,2);
else
    h = linspace(h0, h1, n+is_circ).';
end

% Lose the last hue if it is a whole circumnavigation
if is_circ
    h = h(1:end-1,:);
end

% Make Lab values for all co-ordinates
% c is row vector representing the chroma for each LL
% h is column vector
aa = bsxfun(@times, cc, cosd(h));
bb = bsxfun(@times, cc, sind(h));

% -------------------------------------------------------------------------
% Now pick out these colours

LL_Lab = cell(numel(LL),1);
LL_rgb = cell(numel(LL),1);

for j=1:numel(LL)
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
    I = bsxfun(@plus,(1:n),n*(0:(numel(LL)-1))');
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
    
    % Read the Lab values but leave a NaN between each LL
    LL_Lab2 = nan((n+1)*numel(LL),3);
    I = setdiff(1:((n+1)*numel(LL)),(n+1)*(1:numel(LL)));
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