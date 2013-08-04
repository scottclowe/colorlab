function Lab = gd_rgb2lab(rgb, use_uplab, spacefun)
% Utility for converting from sRGB to CIELab (or UPLab)
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
    Lab = spacefun(rgb);
    return;
end

% If ImageProcessingToolbox is available to use, use it.
% Move from sRGB to CIELab
if license('checkout','image_toolbox')
    Lab = applycform(rgb, makecform('srgb2lab'));
else
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
    % Move from sRGB to CIELab
    Lab = colorspace('RGB->Lab',rgb);
end

% Move from CIELab to UPLab
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
    % A2B0: CIELab to UPLab
    cform = makecform('CLUT', P, 'AToB0');
    Lab = applycform(Lab, cform);
end

end