function rgb = gd_lab2rgb(Lab, use_uplab, spacefun)
% Utility for converting from CIELab (or UPLab) to sRGB
% Gracefully degrading utility
% Uses 'spacefun' if provided ...
% If not, uses ImageProcessingToolbox if available ...
% If not, uses colorspace (which is a function available on FEX) ...
% If not, recommends download of colorspace (using suggestFEXpackage)

if nargin<3
    spacefun = [];
end
if nargin<2
    use_uplab = false;
end

if ~isempty(spacefun)
    rgb = spacefun(Lab);
    return;
end

% Move from UPLab to CIELab
% UPLab was made by Bruce Lindbloom and provides a space where Munsell
% colors are uniformly spaced
% http://www.brucelindbloom.com/index.html?UPLab.html
% The A2B0 tag of the profile converts from CIE Lab to UP Lab. The B2A0 tag 
% performs the inverse transformation... In this way, a linear gamut 
% mapping in UP Lab is equivalent to a nonlinear mapping along the curved 
% Munsell lines in CIE Lab. This should fix, or at least greatly reduce the 
% "blue turns purple" problem.
if use_uplab;
    P = iccread('CIELab_to_UPLab.icc');
    % B2A0: UPLab to CIELab
    cform = makecform('CLUT', P, 'BToA0');
    Lab = applycform(Lab, cform);
end

% If ImageProcessingToolbox is available to use, use it.
% Move from CIELab to sRGB
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
% Move from CIELab to sRGB
rgb = colorspace('Lab->RGB',Lab);
return;

end