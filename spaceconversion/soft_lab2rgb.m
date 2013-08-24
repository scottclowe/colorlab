function rgb = soft_lab2rgb(varargin)
% Utility for converting from CIELab (or UPLab) to sRGB

rgb = hard_lab2rgb(varargin{:});

rgb(rgb<0) = 0;
rgb(rgb>1) = 1;

end