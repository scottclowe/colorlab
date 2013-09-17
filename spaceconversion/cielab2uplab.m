function Lab = cielab2uplab(Lab)

% Move from CIELab to UPLab
% UPLab was made by Bruce Lindbloom and provides a color space where the
% Munsell colors are uniformly spaced.
% http://www.brucelindbloom.com/index.html?UPLab.html
% Turns out, the mapping is the opposite way around to what Bruce says, and
% the original (CIELab_to_UPLab.icc) is and the second version
% (CIELab_to_UPLab2.icc) is not compliant with the standard directions.
% http://argyllcms.com/doc/iccgamutmapping.html

P = iccread('CIELab_to_UPLab.icc');
% BToA0: CIELab to UPLab
Lab = applycform(Lab, makecform('CLUT', P, 'BToA0'));

end