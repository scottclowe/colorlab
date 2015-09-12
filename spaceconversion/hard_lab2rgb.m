function rgb = hard_lab2rgb(Lab, use_uplab)
% Utility for converting from CIELab (or UPLab) to sRGB
% Uses a modified version of colorspace which does not limit the RGB output
% to be between [0,1]

if nargin<2
    use_uplab = false;
end

if use_uplab;
    Lab = uplab2cielab(Lab);
end

rgb = clab_colorspace('Lab->RGB',Lab);

end