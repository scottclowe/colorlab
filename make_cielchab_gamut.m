function [gamut] = make_cielchab_gamut(space, N, point_method, use_uplab)

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
    point_method = 'face-plus';
end
if nargin<4
    use_uplab = false;
end

% - INPUT SWITCHING -
switch lower(space)
    case {'srgb','rgb'}
        space = 'srgb'; % Canonical
    otherwise
        error('Script undefined for space %s. Only sRGB currently supported.',...
            space);
end
switch lower(point_method)
    case {'cube','point'}
        point_method = 'cube'; % Canonical
    case {'face','faces'}
        point_method = 'face'; % Canonical
    case {'face-plus'}
        point_method = 'face-plus'; % Canonical
    otherwise
        error('Unknown point picking method');
end
if isempty(N)
    switch lower(point_method)
        case 'cube'
            N = 256;  % 8-bit filled cube
        case 'face'
            N = 1024; % 64-bit open cube
        case 'face-plus'
            N = 1024;
        otherwise
            error('Unknown point method');
    end
end
% Would use this for CMYK, but not yet set up
% Have to go CMYK->RGB->LAB
% kform = makecform('cmyk2srgb', 'RenderingIntent', 'RelativeColorimetric');

% Make the main gamut object
gamut = make_gamut_lh(space, N, point_method, use_uplab);

% Have to compute again to find a format for the gamut which lets us browse
% by chroma
gamut.lch_chr = make_gamut_chr(gamut, use_uplab);

% Make a mesh we can interpolate on
gamut.lchmesh = make_gamut_mesh(gamut);

end


% Make a gamut indexed by lightness and hue
function gamut = make_gamut_lh(space, N, point_method, use_uplab)

% Don't want these to be 3,6,7,9,11,...
% Only products of 2 and 5 which divide 10 nicely
switch lower(point_method)
    
    case 'cube' % Canonical
        Lintv = 1;
        Lprec = 1;
        hintv = 1;
        hprec = 1;
        
    case 'face' % Canonical
        Lintv = 1/ceil(N/1024);
        Lprec = 1/ceil(N/1024);
        hintv = 1/ceil(N/4096);
        hprec = 1/ceil(N/4096);
        
    case 'face-plus'
        Lintv = 1/ceil(N/1024);
        Lprec = 1/ceil(N/512 );
        hintv = 1/ceil(N/4096);
        hprec = 1/ceil(N/1024);
        
end

switch lower(point_method)
    
    case 'cube' % Canonical
        warning('Cube is inefficient');
        x = linspace(0,1,N);
        [X,Y,Z] = meshgrid(x,x,x);
        rgb = [X(:) Y(:) Z(:)];
        Ntot = N.^3;
        
    case 'face' % Canonical
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
        
    case 'face-plus'
        % Add bonus points for hard-to-reach parts near 0 and 100 lightness
        top_cut1 = 0.04;
        top_cut2 = 0.01;
        ntop = round(N/5);
        nmid = N;
        % Middle
        Nmid = 6*nmid.^2;
        rgbmid = nan(Nmid,3);
        x = linspace(0,1,nmid);
        [X,Y] = meshgrid(x,x);
        rgb0 = [X(:) Y(:) zeros(nmid.^2,1)];
        rgb1 = [X(:) Y(:)  ones(nmid.^2,1)];
        for i=0:2
            j = i*2;
            rgbmid(j*nmid.^2+1:(j+1)*nmid.^2, :) = circshift(rgb0,[0 i]);
            j = i*2+1;
            rgbmid(j*nmid.^2+1:(j+1)*nmid.^2, :) = circshift(rgb1,[0 i]);
        end
        % Top and Bottom 4%
        Ntop = 6*ntop.^2;
        rgbtop1 = nan(Ntop,3);
        x = linspace(0,top_cut1,ntop);
        [X,Y] = meshgrid(x,x);
        rgb0 = [X(:) Y(:) zeros(ntop.^2,1)];
        x = linspace(1-top_cut1,1,ntop);
        [X,Y] = meshgrid(x,x);
        rgb1 = [X(:) Y(:)  ones(ntop.^2,1)];
        for i=0:2
            j = i*2;
            rgbtop1(j*ntop.^2+1:(j+1)*ntop.^2, :) = circshift(rgb0,[0 i]);
            j = i*2+1;
            rgbtop1(j*ntop.^2+1:(j+1)*ntop.^2, :) = circshift(rgb1,[0 i]);
        end
        % Top and Bottom 1%
        Ntop = 6*ntop.^2;
        rgbtop2 = nan(Ntop,3);
        x = linspace(Lprec,top_cut2,ntop);
        [X,Y] = meshgrid(x,x);
        rgb0 = [X(:) Y(:) zeros(ntop.^2,1)];
        x = linspace(1-top_cut2,1-Lprec,ntop);
        [X,Y] = meshgrid(x,x);
        rgb1 = [X(:) Y(:)  ones(ntop.^2,1)];
        for i=0:2
            j = i*2;
            rgbtop2(j*ntop.^2+1:(j+1)*ntop.^2, :) = circshift(rgb0,[0 i]);
            j = i*2+1;
            rgbtop2(j*ntop.^2+1:(j+1)*ntop.^2, :) = circshift(rgb1,[0 i]);
        end
        % Concatenate
        rgb = [rgbmid; rgbtop1; rgbtop2];
        Ntot = size(rgb,1);
        
    otherwise
        error('Unknown point picking method');
end

% Swap between spaces. Use whichever function is available
Lab = gd_rgb2lab(rgb, use_uplab);

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
L = round(L/Lprec)*Lprec;
h = round(h/hprec)*hprec;

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
rgb_gamut = gd_lab2rgb(lab_gamut, use_uplab);

gamut.lch           = lch_gamut;
gamut.lab           = lab_gamut;
gamut.rgb           = rgb_gamut;
gamut.space         = space;
gamut.N             = N;
gamut.Ntot          = Ntot;
gamut.point_method  = point_method;
gamut.Lprec         = Lprec;
gamut.Lintv         = Lintv;
gamut.hprec         = hprec;
gamut.hintv         = hintv;

end


% Make a gamut indexed by chroma
function [lch_chr] = make_gamut_chr(g, use_uplab)

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
Lmin.rgb = gd_lab2rgb(Lmin.lab, use_uplab);
Lmax.rgb = gd_lab2rgb(Lmax.lab, use_uplab);


lch_chr.Lmin = Lmin;
lch_chr.Lmax = Lmax;

end
