function cmap = cie_rainbow_cmap_make(n, attributes, spacefun, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<4 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end
if nargin<3
    spacefun = []; % function to map from cielab to srgb
end
if nargin<2 || isempty(attributes)
    attributes = 'greenmid'; % Colormap type option
end
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1); % Number of colours in the colormap
end

% -------------------------------------------------------------------------
% Main function body
attrbreaks = strfind(attributes,':');
if isempty(attrbreaks)
    typ = attributes;
    paramset = '1';
else
    typ = attributes(1:attrbreaks(1)-1);
    paramset = attributes(attrbreaks(1)+1:end);
end
    
switch lower(typ)
    case 'smooth'
        cmap = cie_rainbow_cmap_make_smooth(n, paramset, spacefun, dbg);
    case 'greenmid'
        cmap = cie_rainbow_cmap_make_greenmid(n, paramset, spacefun, dbg);
    otherwise
        error('Unfamiliar colormap attribute: %s',typ);
end

end


function params = get_rainbow_ellipse_params(paramset)
% Declares parameters for the (partial) ellipse in CIELab which is used to
% generate the colormap

switch paramset
    case '1'
        params.use_uplab = false;
        params.P       = [-60 49 88];                       % Green point
        params.C       = [301.0341 -123.2044  -32.2823];    % Centre
        params.a       = 400;                               % Major axis coef
        params.b       = 106;                               % Minor axis coef
        params.theta   = -25.5/180*pi; % Rotation of ellipse wrt ab-plane of L*a*b*
        params.phi     = -17.5/180*pi; % L* Inclination of major axis
        params.psi     =   1.7/180*pi; % L* Inclination of minor axis
        params.U       = [ cos(params.theta) sin(params.theta) sin(params.phi)]; % Major axis
        params.V       = [-sin(params.theta) cos(params.theta) sin(params.psi)]; % Minor axis
        params.t_start = 4.2;          % Good start point
        params.t_end   = 2.3;          % Good end point

    case '2'
        params.use_uplab = false;
        params.P       = [-41 37 91];
        params.C       = [131.8926  -41.7917   15.2377];
        params.a       = 190;
        params.b       =  89;
        params.theta   = -24.5/180*pi;
        params.phi     = -23.5/180*pi;
        params.psi     =     0/180*pi;
        params.U       = [ cos(params.theta) sin(params.theta) sin(params.phi)]; % Major axis
        params.V       = [-sin(params.theta) cos(params.theta) sin(params.psi)]; % Minor axis
        params.t_start = 4.6;          % Good start point
        params.t_end   = 2.05;         % Good end point
        
    otherwise
        error('Unfamiliar parameter set');
end

end


function [btsp_t_srt, btsp_t_end] = find_ellipse_ends(paramset)
% Finds the limiting points which are only just inside the ellipse

params = get_rainbow_ellipse_params(paramset);
rgbgamut = fetch_cielchab_gamut('srgb', 2048, 'face');


% First generate a ton of points to find the end points of theta
n = 100001;
t = linspace(params.t_start, params.t_end, n);
x = params.C(1) + params.a * cos(t) * params.U(1) + params.b * sin(t) * params.V(1);
y = params.C(2) + params.a * cos(t) * params.U(2) + params.b * sin(t) * params.V(2);
z = params.C(3) + params.a * cos(t) * params.U(3) + params.b * sin(t) * params.V(3);

Lab = [z' x' y'];
TF = isingamut(Lab,rgbgamut,'Lab');

% Find edges of gamut
I_srt = find(TF,1,'first');
I_end = find(TF,1,'last');

% Go back by a Euclidian distance of 1 to be sure we're in gamut
de_srt = sqrt(sum(bsxfun(@minus,Lab(I_srt:I_end,:),Lab(I_srt,:)).^2,2));
de_end = sqrt(sum(bsxfun(@minus,Lab(I_srt:I_end,:),Lab(I_end,:)).^2,2));
btsp_t_srt = t(I_srt-1+find(de_srt>1,1,'first'));
btsp_t_end = t(I_srt-1+find(de_end>1,1,'last'));

end


function ciebow_cmap = cie_rainbow_cmap_make_smooth(n_target, paramset, spacefun, dbg)

params = get_rainbow_ellipse_params(paramset);
rgbgamut = fetch_cielchab_gamut('srgb', 2048, 'face');

% First generate a ton of points to find the end points of theta
[btsp_t_srt, btsp_t_end] = find_ellipse_ends(paramset);

% Now generate a ton of points in this range
n = 100001;
t = linspace(btsp_t_srt, btsp_t_end, n);
x = params.C(1) + params.a * cos(t) * params.U(1) + params.b * sin(t) * params.V(1);
y = params.C(2) + params.a * cos(t) * params.U(2) + params.b * sin(t) * params.V(2);
z = params.C(3) + params.a * cos(t) * params.U(3) + params.b * sin(t) * params.V(3);

Lab = [z' x' y'];
TF = isingamut(Lab,rgbgamut,'Lab');

% Check all points are now inside
% if dbg;
%     fprintf('This need to be equal %d = %d\t%d\n',sum(TF),length(TF),(sum(TF)/length(TF)));
% end
if sum(TF)~=length(TF)
    warning('Not all points in gamut');
end

% Using simple 1931 color difference metric
% Find Euclidian 2norm distance between each set of consequtive points
Lab_dif = diff(Lab,1,1);
Lab_inv = sqrt(sum(Lab_dif.^2,2));

% % % Using more up to date color difference metric
% % % https://en.wikipedia.org/wiki/Color_difference#CIEDE2000
% %
% % My implementation (doesn't work)
% % LCh = applycform(Lab,makecform('lab2lch'));
% % Lab_inv = ciede(LCh(1:end-1,:),LCh(2:end,:));
% %
% % From Westland and Ripamonti (2004) - FEX:40640
% Lab_inv = cie00de(Lab(1:end-1,:),Lab(2:end,:));

% Hence find length of complete arc
arclen_full = sum(Lab_inv);

% Find 256 points equally spaced out of the 100000 we generated
% n_target = 256;
arclen_inv_target = arclen_full./(n_target-1);

% Match target to the arclength NOT the euclidean length between points
% This is easier to code
Lab_inv_runtot = [0; cumsum(Lab_inv)];
arclen_targets = arclen_inv_target * (0:(n_target-1));

t_new = interp1(Lab_inv_runtot, t, arclen_targets);

t_new(1)   = btsp_t_srt;
t_new(end) = btsp_t_end;

x = params.C(1) + params.a * cos(t_new) * params.U(1) + params.b * sin(t_new) * params.V(1);
y = params.C(2) + params.a * cos(t_new) * params.U(2) + params.b * sin(t_new) * params.V(2);
z = params.C(3) + params.a * cos(t_new) * params.U(3) + params.b * sin(t_new) * params.V(3);

Lab = [z' x' y'];
ciebow_cmap = gd_lab2rgb(Lab, params.use_uplab, spacefun);

% debugging figures
if dbg;
    img = repmat(ciebow_cmap,[1 1 20]);
    img = permute(img,[1 3 2]);

    figure;
    imagesc(img);
    axis xy;

    [TF,P2] = isingamut(Lab,rgbgamut,'Lab');

    % disp(sum(TF)/length(TF))

    figure;
    hold on;
    plot3(Lab(:,2), Lab(:,3), Lab(:,1), 'bx')
    plot3(P2(:,2) , P2(:,3) , P2(:,1) , 'r-')
    xlabel('a');
    ylabel('b');
    zlabel('L');
    view(0,90);
    title('Rotate to see L');
    
    
    figure;
    hold on;
    plot3(Lab(:,2), Lab(:,3), Lab(:,1), 'kx-');
    
    % Get a mesh version of the gamut
    if ~isfield(rgbgamut,'lchmesh')
        rgbgamut.lchmesh = make_gamut_mesh(rgbgamut);
    end

    L = rgbgamut.lchmesh.Lgrid([1:4:(end-1) 1],[1:4:(end-1) end]);
    c = rgbgamut.lchmesh.cgrid([1:4:(end-1) 1],[1:4:(end-1) end]);
    h = rgbgamut.lchmesh.hgrid([1:4:(end-1) 1],[1:4:(end-1) end])/180*pi;
    a = c.*cos(h);
    b = c.*sin(h);
    
    CMAP = gd_lab2rgb([L(:) a(:) b(:)], use_uplab);

    hs = mesh(a,b,L,reshape(CMAP,[size(L) 3]));
    set(hs,'FaceColor','none');

    set(gca,'Color',[0.4663 0.4663 0.4663]);
    set(gca,'XLim',[-150 150],'YLim',[-150 150],'ZLim',[0 100]);
    xlabel('a*')
    ylabel('b*')
    zlabel('L*')
    
    
    plot_labcurve_rgbgamut(Lab)
end

end


function ciebow_cmap = cie_rainbow_cmap_make_greenmid(n_target, paramset, spacefun, dbg)
% Option to have same amount of red as blue
% This is a very simple method used where the curve is split in two and
% there are different Delta E values between colours in each.
% It would be better to have the rate change sigmoidally instead

% Mid point is at t=pi

params = get_rainbow_ellipse_params(paramset);
rgbgamut = fetch_cielchab_gamut('srgb', 2048, 'face');

% First generate a ton of points to find the end points of theta
[btsp_t_srt, btsp_t_end] = find_ellipse_ends(paramset);

% Now generate a ton of points in this range
n = 100001;
t1 = linspace(btsp_t_srt, pi, n);
x1 = params.C(1) + params.a * cos(t1) * params.U(1) + params.b * sin(t1) * params.V(1);
y1 = params.C(2) + params.a * cos(t1) * params.U(2) + params.b * sin(t1) * params.V(2);
z1 = params.C(3) + params.a * cos(t1) * params.U(3) + params.b * sin(t1) * params.V(3);
t2 = linspace(pi, btsp_t_end, n);
x2 = params.C(1) + params.a * cos(t2) * params.U(1) + params.b * sin(t2) * params.V(1);
y2 = params.C(2) + params.a * cos(t2) * params.U(2) + params.b * sin(t2) * params.V(2);
z2 = params.C(3) + params.a * cos(t2) * params.U(3) + params.b * sin(t2) * params.V(3);

Lab1 = [z1' x1' y1'];
Lab2 = [z2' x2' y2'];


% Check all points are now inside
TF = isingamut([Lab1;Lab2],rgbgamut,'Lab');
% if dbg;
%     fprintf('This need to be equal %d = %d\t%d\n',sum(TF),length(TF),(sum(TF)/length(TF)));
% end
if sum(TF)~=length(TF)
    warning('Not all points in gamut');
end

% Using simple 1931 color difference metric
% Find Euclidian 2norm distance between each set of consequtive points
Lab1_dif = diff(Lab1,1,1);
Lab2_dif = diff(Lab2,1,1);
Lab1_inv = sqrt(sum(Lab1_dif.^2,2));
Lab2_inv = sqrt(sum(Lab2_dif.^2,2));

% % % Using more up to date color difference metric
% % % https://en.wikipedia.org/wiki/Color_difference#CIEDE2000
% %
% % My implementation (doesn't work)
% % LCh1 = applycform(Lab1,makecform('lab2lch'));
% % LCh2 = applycform(Lab2,makecform('lab2lch'));
% % Lab1_inv = ciede(LCh1(1:end-1,:),LCh1(2:end,:)); % Coded by me
% % Lab2_inv = ciede(LCh2(1:end-1,:),LCh2(2:end,:)); % Coded by me
% %
% % From Westland and Ripamonti (2004) - FEX:40640
% Lab1_inv = cie00de(Lab1(1:end-1,:),Lab1(2:end,:)); % From Westland and Ripamonti (2004)
% Lab2_inv = cie00de(Lab2(1:end-1,:),Lab2(2:end,:)); % From Westland and Ripamonti (2004)


% Hence find length of complete arc
arclen_blu = sum(Lab1_inv);
arclen_red = sum(Lab2_inv);

% Find 256 points equally spaced out of the 100000 we generated
% n_target = 256;
if mod(n_target,2);
    n_each = (n_target-1)/2;
else
    n_each = n_target/2;
end
arclen_inv_target_blu = arclen_blu./n_each;
arclen_inv_target_red = arclen_red./n_each;

% Match target to the arclength NOT the euclidean length between points
% This is easier to code
Lab1_inv_runtot = [0; cumsum(Lab1_inv)];
Lab2_inv_runtot = [0; cumsum(Lab2_inv)];
arclen_targets_blu = arclen_inv_target_blu * (0:n_each);
arclen_targets_red = arclen_inv_target_red * (0:n_each);

t1_new = interp1(Lab1_inv_runtot, t1, arclen_targets_blu);
t2_new = interp1(Lab2_inv_runtot, t2, arclen_targets_red);

% If even number requested, there is an extra colour. Sorry.
t_new = [btsp_t_srt t1_new(2:end-1) pi t2_new(2:end-1) btsp_t_end];

x = params.C(1) + params.a * cos(t_new) * params.U(1) + params.b * sin(t_new) * params.V(1);
y = params.C(2) + params.a * cos(t_new) * params.U(2) + params.b * sin(t_new) * params.V(2);
z = params.C(3) + params.a * cos(t_new) * params.U(3) + params.b * sin(t_new) * params.V(3);

Lab = [z' x' y'];
ciebow_cmap = gd_lab2rgb(Lab, params.use_uplab, spacefun);

% debugging figures
if dbg;
    img = repmat(ciebow_cmap,[1 1 20]);
    img = permute(img,[1 3 2]);

    figure;
    imagesc(img);
    axis xy;

    [TF,P2] = isingamut(Lab,rgbgamut,'Lab');

    % disp(sum(TF)/length(TF))

    figure;
    hold on;
    plot3(Lab(:,2), Lab(:,3), Lab(:,1), 'kx')
    plot3(P2(:,2) , P2(:,3) , P2(:,1) , 'r-')
    xlabel('a');
    ylabel('b');
    zlabel('L');
    view(0,90);
    title('Rotate to see L')
    
    
    figure;
    hold on;
    plot3(Lab(:,2), Lab(:,3), Lab(:,1), 'kx-');
    
    % Get a mesh version of the gamut
    if ~isfield(rgbgamut,'lchmesh')
        rgbgamut.lchmesh = make_gamut_mesh(rgbgamut);
    end

    L = rgbgamut.lchmesh.Lgrid([1:4:(end-1) 1],[1:4:(end-1) end]);
    c = rgbgamut.lchmesh.cgrid([1:4:(end-1) 1],[1:4:(end-1) end]);
    h = rgbgamut.lchmesh.hgrid([1:4:(end-1) 1],[1:4:(end-1) end])/180*pi;
    a = c.*cos(h);
    b = c.*sin(h);

    CMAP = gd_lab2rgb([L(:) a(:) b(:)], use_uplab);

    hs = mesh(a,b,L,reshape(CMAP,[size(L) 3]));
    set(hs,'FaceColor','none');

    set(gca,'Color',[0.4663 0.4663 0.4663]);
    set(gca,'XLim',[-150 150],'YLim',[-150 150],'ZLim',[0 100]);
    xlabel('a*')
    ylabel('b*')
    zlabel('L*')
    
    plot_labcurve_rgbgamut(Lab)
end

end