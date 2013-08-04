function rgb = gd_lab2rgb(Lab, func)
% Utility for converting from CIELab to sRGB
% Gracefully degrading utility
% Uses 'func' if provided ...
% If not, uses ImageProcessingToolbox if available ...
% If not, uses colorspace (which is a function available on FEX) ...
% If not, recommends download of colorspace (using suggestFEXpackage)

if nargin<2
    func = [];
end

if ~isempty(func)
    rgb = func(Lab);
    return;
end

% If ImageProcessingToolbox is available to use, use it.
if license('checkout','image_toolbox')
    rgb = applycform(Lab, makecform('lab2srgb'));
    return;
end

% If we don't have colorspace, suggest it
if ~exist('colorspace','file')
    if exist('suggestFEXpackage','file')
        folder = suggestFEXpackage(28790,...
            ['Since the Image Processing Toolbox is unavailable, '...
             'you may wish to download the colorspace package.\n' ...
             'This package will allow you to convert between different '...
             'colorspaces without the MATLAB toolbox' ...
            ]);
    else
        folder = '';
    end
    if isempty(folder)
        % User did not download colorspace
        error('LABWHEEL:NoIPToolbox:NoColorspace',...
            ['Could not checkout the Image Processing Toolbox. ' ...
             'Colorspace function not present either.\n' ...
             'You need one of the two to run this function. ' ...
             ]);
    end
end

% Use colorspace
rgb = colorspace('Lab->RGB',Lab);
return;

end