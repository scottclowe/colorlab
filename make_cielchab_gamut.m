function [gamut] = make_cielchab_gamut(space, N, point_method)

% - INPUT HANDLING -
if nargin<1
    space = 'srgb';
end
if nargin<2
    % Number of data points will be N^3   in cube  method
    % Number of data points will be 6*N^2 in faces method
    N = []; 
end
if nargin<3
    point_method = 'face';
end
if isempty(N)
    switch lower(point_method)
        case {'cube','point'}
            N = 256;  % 8-bit filled cube
        case {'face'}
            N = 1024; % 64-bit open cube
        otherwise
            error('Unknown point method');
    end
end

% - INPUT SWITCHING -
Lprec = ceil(N/1024);
Lintv = 1/Lprec;
Lprec = Lprec/2;
hprec = ceil(N/2048);
hintv = 1;
switch lower(space)
    case {'srgb','rgb'}
        space = 'srgb'; % Canonical
    otherwise
        error('Script undefined for space %s. Only sRGB currently supported.',...
            space);
end

% Would use this for CMYK, but not yet set up
% Have to go CMYK->RGB->LAB
% kform = makecform('cmyk2srgb', 'RenderingIntent', 'RelativeColorimetric');

switch lower(point_method)
    case {'cube','point'}
        point_method = 'cube'; % Canonical
        x = linspace(0,1,N);
        [X,Y,Z] = meshgrid(x,x,x);
        rgb = [X(:) Y(:) Z(:)];
        Ntot = N.^3;
    case {'face','faces'}
        point_method = 'face'; % Canonical
        Ntot = 6*N.^2;
        rgb = nan(Ntot,3);
        x = linspace(0,1,N);
        [X,Y] = meshgrid(x,x);
        rgb0 = [X(:) Y(:) zeros(N.^2,1)];
        rgb1 = [X(:) Y(:)  ones(N.^2,1)];
        for i=0:2
            j = i*2;
            rgb(j*N.^2+1:(j+1)*N.^2, :) = circshift(rgb0,[0 i]);
            j = i*2+1;
            rgb(j*N.^2+1:(j+1)*N.^2, :) = circshift(rgb1,[0 i]);
        end
    otherwise
        error('Unknown point picking method');
end

% Swap between spaces. Use whichever function is available
if license('checkout','image_toolbox')
    % If using ImageProcessingToolbox
    cform = makecform('srgb2lab');
    Lab = applycform(rgb, cform);
elseif exist('colorspace','file')
    % Use colorspace
%     warning('LABWHEEL:NoIPToolbox:UseColorspace',...
%         ['Could not checkout the Image Processing Toolbox. ' ...
%          'Using colorspace function, but output not guaranteed correct.']);
    Lab = colorspace('RGB->Lab',rgb);
else
    % Use colorspace
    warning('LABWHEEL:NoIPToolbox:NoColorspace',...
        ['Could not checkout the Image Processing Toolbox. ' ...
         'Colorspace function not present either.\n' ...
         'You need one of the two to run this function.']);
     if exist('suggestFEXpackage','file')
        suggestFEXpackage(28790,'You may wish to download the colorspace package.');
     end
end

% Seperate Lab components
L = Lab(:,1);
a = Lab(:,2);
b = Lab(:,3);

% Sort in order of ascending lightness
[L,I] = sort(L,1,'ascend');
a = a(I);
b = b(I);

% Find CIE LCh_ab
% chroma is radial distance from cylindrical pole (pure grey line)
c = (a.^2 + b.^2).^(1/2);
% hue is angle anticlockwise from the a* axis
% atan has a period of pi, so we need to do some jigging around to preserve
% accurate angle
% h = atan(b./a); % In radians
% h(isnan(h)) = 0;
% h = h./(2*pi)*360; % In degrees
% s = sign(a);
% s(s==1) = 0;
% h = h-180*s;
% h = mod(h,360);
% Should use atan2 instead!
h = mod(atan2(b,a)/pi*180,360);

% Round to nearest (integer/Lprec) for lightness and to nearest degree for hue
L = round(L*Lprec)/Lprec;
h = round(h*hprec)/hprec;

% Slice the space into layers of different lightness
% For each slice, find the maximum chroma at each angle of hue
% Add each datapoint to the list defining the edges of the gamut
LLL = 0:Lintv:100;
hhh = 0:hintv:359.999;
lch_gamut = nan(length(LLL)*length(hhh),3);
% Set all hues for L=0 manually
lch_gamut(1:length(hhh),[1 2]) = zeros(length(hhh),2);
lch_gamut(1:length(hhh),3)     = hhh;
% Now process the middle
i = length(hhh);
for LL = setdiff(LLL,[0 100])
    vi = L==LL;
    LLh = h(vi);
    LLc = c(vi);
    for hh = hhh
        % Find maximum chroma at this Lightness and Hue
        cc = max(LLc(LLh==hh));
        if isempty(cc); continue; end
        i = i+1;
        lch_gamut(i,:) = [LL cc hh];
    end
end
% Set all hues for L=100 manually
lch_gamut(i+(1:length(hhh)),1) = 100*ones(length(hhh),1);
lch_gamut(i+(1:length(hhh)),2) = zeros(length(hhh),1);
lch_gamut(i+(1:length(hhh)),3) = hhh;
i = i+length(hhh);
% Chop off the excess points which we couldn't find values for
lch_gamut = lch_gamut(1:i,:);

ga = lch_gamut(:,2).*cos(lch_gamut(:,3)/360*(2*pi));
gb = lch_gamut(:,2).*sin(lch_gamut(:,3)/360*(2*pi));
lab_gamut = [lch_gamut(:,1) ga gb];

% Move back to RGB so we have a set of colors we can show
if license('checkout','image_toolbox')
    % If using ImageProcessingToolbox
    cform = makecform('lab2srgb');
    rgb_gamut = applycform(lab_gamut, cform);
elseif exist('colorspace','file')
    % Use colorspace
%     warning('LABWHEEL:NoIPToolbox:UseColorspace',...
%         ['Could not checkout the Image Processing Toolbox. ' ...
%          'Using colorspace function, but output not guaranteed correct.']);
    rgb_gamut = colorspace('Lab->RGB',lab_gamut);
else
    % Use colorspace
    warning('LABWHEEL:NoIPToolbox:NoColorspace',...
        ['Could not checkout the Image Processing Toolbox. ' ...
         'Colorspace function not present either.\n' ...
         'You need one of the two to run this function.']);
     suggestFEXpackage(28790,'You may wish to download the colorspace package.')
end

gamut.lch           = lch_gamut;
gamut.lab           = lab_gamut;
gamut.rgb           = rgb_gamut;
gamut.space         = space;
gamut.N             = N;
gamut.point_method  = point_method;
gamut.Lprec         = Lprec;
gamut.Lintv         = Lintv;
gamut.hprec         = hprec;
gamut.hintv         = hintv;


% Have to compute again to find a format for the gamut which lets us browse
% by chroma
gamut.lch_chr = find_gamut_chr(gamut);

end


function [lch_chr] = find_gamut_chr(g)
max_c = max(g.lch(:,3));
hues = unique(g.lch(:,3));
Lmin.lch = nan((max_c+1)*length(hues),3);
Lmax.lch = nan((max_c+1)*length(hues),3);
i = 0;
for c=0:max_c
    for ihue = 1:length(hues)
        h = hues(ihue);
        sublch = g.lch(g.lch(:,3)==h,:);
        sublch = sublch(sublch(:,2)>=c,:);
        if isempty(sublch); continue; end
        i = i+1;
        [C,I] = min(sublch(:,1));
        Lmin.lch(i,:) = [sublch(I,1) c h];
        [C,I] = max(sublch(:,1));
        Lmax.lch(i,:) = [sublch(I,1) c h];
    end
end
Lmin.lch = Lmin.lch(1:i,:);
Lmax.lch = Lmax.lch(1:i,:);

ga = Lmin.lch(:,2).*cos(Lmin.lch(:,3)/360*(2*pi));
gb = Lmin.lch(:,2).*sin(Lmin.lch(:,3)/360*(2*pi));
Lmin.lab = [Lmin.lch(:,1) ga gb];

ga = Lmax.lch(:,2).*cos(Lmax.lch(:,3)/360*(2*pi));
gb = Lmax.lch(:,2).*sin(Lmax.lch(:,3)/360*(2*pi));
Lmax.lab = [Lmax.lch(:,1) ga gb];

% Move back to RGB so we have a set of colors we can show
if license('checkout','image_toolbox')
    % If using ImageProcessingToolbox
    cform = makecform('lab2srgb');
    Lmin.rgb = applycform(Lmin.lab, cform);
    Lmax.rgb = applycform(Lmax.lab, cform);
elseif exist('colorspace','file')
    % Use colorspace
%     warning('LABWHEEL:NoIPToolbox:UseColorspace',...
%         ['Could not checkout the Image Processing Toolbox. ' ...
%          'Using colorspace function, but output not guaranteed correct.']);
    Lmin.rgb = colorspace('Lab->RGB',Lmin.lab);
    Lmax.rgb = colorspace('Lab->RGB',Lmax.lab);
else
    % Use colorspace
    warning('LABWHEEL:NoIPToolbox:NoColorspace',...
        ['Could not checkout the Image Processing Toolbox. ' ...
         'Colorspace function not present either.\n' ...
         'You need one of the two to run this function.']);
     suggestFEXpackage(28790,'You may wish to download the colorspace package.')
end

lch_chr.Lmin = Lmin;
lch_chr.Lmax = Lmax;

end
