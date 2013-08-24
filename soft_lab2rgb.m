function rgb = soft_lab2rgb(Lab, use_uplab)
% Utility for converting from CIELab (or UPLab) to sRGB
% Gracefully degrading utility
% Uses 'spacefun' if provided ...
% If not, uses ImageProcessingToolbox if available ...
% If not, uses colorspace (which is a function available on FEX) ...
% If not, recommends download of colorspace (using suggestFEXpackage)

if nargin<2
    use_uplab = false;
end


% Move from UPLab to CIELab
% UPLab was made by Bruce Lindbloom and provides a color space where the
% Munsell colors are uniformly spaced.
% http://www.brucelindbloom.com/index.html?UPLab.html
% Turns out, the mapping is the opposite way around to what Bruce says, and
% the original (CIELab_to_UPLab.icc) is and the second version
% (CIELab_to_UPLab2.icc) is not compliant with the standard directions.
% http://argyllcms.com/doc/iccgamutmapping.html
if use_uplab;
    P = iccread('CIELab_to_UPLab.icc');
    % AToB0: UPLab to CIELab
    cform = makecform('CLUT', P, 'AToB0');
    Lab = applycform(Lab, cform);
end

% % If ImageProcessingToolbox is available to use, use it.
% % Move from CIELab to sRGB
% if license('checkout','image_toolbox')
%     rgb = applycform(Lab, makecform('lab2srgb'));
%     return;
% end

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

% rgb = my_colorspace('Lab->RGB',Lab);

end