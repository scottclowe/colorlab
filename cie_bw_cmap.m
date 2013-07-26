function cmap = cie_bw_cmap(n, func)

if nargin<1
    n = size(get(gcf,'colormap'),1);
end
if nargin<2
    func = [];
end

L = linspace(0,100,n)';
a = zeros(n,1);
b = zeros(n,1);

Lab = [L a b];

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
         'Using colorspace function.']);
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