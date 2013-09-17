function Lab = rgb2lab(rgb, use_uplab)
% Utility for converting from sRGB to CIELab (or UPLab)

if nargin<2
    use_uplab = false;
end

Lab = clab_colorspace('RGB->Lab',rgb);

if use_uplab;
    Lab = cielab2uplab(Lab);
end

end