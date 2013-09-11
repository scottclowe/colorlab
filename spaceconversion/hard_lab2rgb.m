function rgb = hard_lab2rgb(Lab, use_uplab)
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
    Lab = applycform(Lab, makecform('CLUT', P, 'AToB0'));
%     tempLab = nan(size(Lab));
%     li = ~(isnan(Lab(:,1)) | isnan(Lab(:,2)) | isnan(Lab(:,3)));
%     tempLab(li,:) = applycform(Lab(li,:), cform);
%     Lab = tempLab;
end

rgb = clab_colorspace('Lab->RGB',Lab);

end