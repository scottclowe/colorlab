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

if use_uplab;
    Lab = uplab2cielab(Lab);
end

rgb = clab_colorspace('Lab->RGB',Lab);

end