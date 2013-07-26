function cmap = cie_blueyellow_cmap(n, func)

if nargin<1
    n = size(get(gcf,'colormap'),1);
end
if nargin<2
    func = [];
end

% CIE       [  L*   a*   b*]

% rgbblue   = [  32   79 -108];
% rgbyellow = [  97  -22   94];

% ciestart = rgbblue;
% cieend = rgbyellow;

ciestart   = [  58   0   -65];
cieend     = [  89   0    87];


Ls = linspace(ciestart(1),cieend(1),n);
as = linspace(ciestart(2),cieend(2),n);
bs = linspace(ciestart(3),cieend(3),n);

Lab = [Ls' as' bs'];

if ~isempty(func)
    cmap = func(Lab);
elseif license('checkout','image_toolbox')
    % If using ImageProcessingToolbox
    cform = makecform('lab2srgb');
    cmap = applycform(Lab, cform);
elseif exist('colorspace','file')
    % Use colorspace
    warning('LABWHEEL:NoIPToolbox:UseColorspace',...
        ['Could not checkout the Image Processing Toolbox. ' ...
         'Using colorspace function, but output not guaranteed correct.']);
    cmap = colorspace('Lab->RGB',Lab);
else
    % Use colorspace
    warning('LABWHEEL:NoIPToolbox:NoColorspace',...
        ['Could not checkout the Image Processing Toolbox. ' ...
         'Colorspace function not present either.\n' ...
         'You need one of the two to run this function.']);
     if exist('suggestFEXpackage','file')
        suggestFEXpackage(28790,'You may wish to download the colorspace package.');
     end
end