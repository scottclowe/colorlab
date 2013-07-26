function Y = find_gamut_xy_maxY(x, y, primaries)
% Could add an option to specify xyY primaries
% Would have to compute rgb2xyz from these
% This is done by converting xyY to XYZ
% Then invert matrix to find xyz2rgb

% primaries should be of the form
% P = [ Red_x Green_x Blue_x;
%       Red_y Green_y Blue_y;
%       Red_Y Green_Y Blue_Y]

% sRGB given below
% These matrices pre-multiply an RGB or XYZ column vector

% http://en.wikipedia.org/wiki/SRGB#Specification_of_the_transformation
% srgb2xyz = [ ...
%     0.4124  0.3576  0.1805;
%     0.2126  0.7152  0.0722;
%     0.0193  0.1192  0.9505];
%
% xyz2srgb = [ ...
%      3.2406 -1.5392 -0.4986;
%     -0.9689  1.8758  0.0415;
%      0.0557 -0.2040  1.0570];

% http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
srgb2xyz = [ ...
 0.4124564  0.3575761  0.1804375;
 0.2126729  0.7151522  0.0721750;
 0.0193339  0.1191920  0.9503041];

xyz2srgb = [ ...
 3.2404542 -1.5371385 -0.4985314;
-0.9692660  1.8760108  0.0415560;
 0.0556434 -0.2040259  1.0572252];

%% Input handling

if nargin<3
    primaries = 'srgb';
end

if ~isequal(size(x),size(y))
    error('Input sizes do not match');
end
shp = size(x);
x = x(:)';
y = y(:)';

if ischar(primaries)
    switch lower(primaries)
        case 'srgb'
            rgb2xyz = srgb2xyz;
%             xyz2rgb = xyz2srgb;
        otherwise
            error('Unknown color space');
    end
else
    rgb2xyz = zeros(size(primaries));
    rgb2xyz(1,:) = primaries(1,:).*primaries(3,:)./primaries(2,:);
    rgb2xyz(2,:) = primaries(3,:);
    rgb2xyz(3,:) = (1-primaries(1,:)-primaries(2,:)).*primaries(3,:)./primaries(2,:);
    
%     xyz2rgb = inv(rgb2xyz);
end

%% Do computation

z = 1-x-y;
% XYZ is scaled version of xyz
% rgb = xyz2rgb * [x;y;z;];
rgb = rgb2xyz \ [x;y;z;]; % Avoid finding matrix inverse
% Find scaling factor which will maximise one of RGB components
Y_over_y = 1./max(rgb, [], 1);
Y = Y_over_y .* y;
% Make anything out of gamut have Y=NaN
Y(any(rgb<0,1)) = NaN;

%% Output handling

Y = reshape(Y,shp);

end