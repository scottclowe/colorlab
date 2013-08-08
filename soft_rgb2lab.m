function Lab = soft_rgb2lab(rgb, use_uplab, spacefun)
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


% Lab = my_colorspace('RGB->Lab',rgb);

% If ImageProcessingToolbox is available to use, use it.
% Move from sRGB to CIELab
if false % license('checkout','image_toolbox')
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
% UPLab was made by Bruce Lindbloom and provides a color space where the
% Munsell colors are uniformly spaced.
% http://www.brucelindbloom.com/index.html?UPLab.html
% Turns out, the mapping is the opposite way around to what Bruce says, and
% the original (CIELab_to_UPLab.icc) is and the second version
% (CIELab_to_UPLab2.icc) is not compliant with the standard directions.
% http://argyllcms.com/doc/iccgamutmapping.html
if use_uplab;
    P = iccread('CIELab_to_UPLab.icc');
    % BToA0: CIELab to UPLab
    cform = makecform('CLUT', P, 'BToA0');
    Lab = applycform(Lab, cform);
end

end