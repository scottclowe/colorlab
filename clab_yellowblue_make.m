function cmap = clab_yellowblue_make(n, attr, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<3 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end
if nargin<2 || isempty(attr)
    attr = ''; % Colormap type option
end
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1); % Number of colours in the colormap
end

% -------------------------------------------------------------------------
use_uplab = false;

switch attr
    case ''
        % CIE  [  L*     C     h ]
%         lch1 = [59.375 63.906  98];
%         lch2 = [59.375 63.906 278];
        
%         lch1 = [64.375 68.8125 102.75];
%         lch2 = [64.375 68.8125 307   ];
        
%         lch1 = [59.125 64.375 103];
%         lch2 = [59.125 64.375 278];
        
        lch1 = [59.5   64.5   102];
        lch2 = [59.5   64.5   282];

    case {'a0','a0:lfix','a0:lfix:1'}
        % CIE  [  L*     C     h ]
        lch1 = [63.375 56.968  90];
        lch2 = [63.375 56.968 270];
        
% OLD
%     case 'a0'
%         % CIE  [  L*    a*  b*]
%         lab1 = [  58    0  -65];
%         lab2 = [  89    0   87];
%         
%     case {'a0:lfix','a0:lfix:1'}
%         % NB: The local maxima for joint-chroma is at h=88 (L=59.5, C=63)
%         % So this is pretty much optimal for fixed L
%         % It used to look like a local maxima, but it sure doesn't anymore
%         % CIE  [  L*    a*  b*]
%         lab1 = [  59.5  0  -62];
%         lab2 = [  59.5  0   62];

% OLDER
%     case 'a0:lfix:2'
%         % CIE  [  L*    a*  b*]
%         lab1 = [ 57.75  0  -65];
%         lab2 = [ 57.75  0   61];
%
%     case 'maxsep'
% %         lch1 = [ 30.5   130    302];
% %         lch2 = [ 90     100    122];
%         lab1 = [ 30.50  68 -109];
%         lab2 = [ 90.00 -52   84];

    case 'primaries'
        
%         lab1 = [97.1279  -21.5484   94.4728]; % yellow
%         lab2 = [32.3006   79.1839 -107.8569]; % blue
        
        lch1 = [97.1279   96.8991  102.8489]; % yellow
        lch2 = [32.3006  133.8029  306.2845]; % blue
        
    otherwise
        error('Unfamiliar colormap attribute: %s',attr);
end

% -------------------------------------------------------------------------
lab1 = [lch1(1) lch1(2).*cosd(lch1(3)) lch1(2).*sind(lch1(3))];
lab2 = [lch2(1) lch2(2).*cosd(lch2(3)) lch2(2).*sind(lch2(3))];

if mod(lch1(3)-lch2(3),360)==180

    L = linspace(lab1(1), lab2(1), n);
    a = linspace(lab1(2), lab2(2), n);
    b = linspace(lab1(3), lab2(3), n);

    Lab = [L' a' b'];

else
    L_grey = (lch1(1)*lch2(2)+lch2(1)*lch1(2))/(lch1(2)+lch2(2));

    neach = floor(n/2)+1;
    % => if n=odd  then neach+(neach-1)=n
    %    if n=even then neach+(neach-1)=n+1

    L1 = linspace(lab1(1), L_grey, neach);
    a1 = linspace(lab1(2),      0, neach);
    b1 = linspace(lab1(3),      0, neach);

    L2 = linspace(L_grey, lab2(1), neach);
    a2 = linspace(     0, lab2(2), neach);
    b2 = linspace(     0, lab2(3), neach);

    Lab1 = [L1' a1' b1'];
    Lab2 = [L2' a2' b2'];
    
    % Convert from Lab to srgb
    Lab  = [Lab1(1:end-1,:); Lab2];
end

cmap = hard_lab2rgb(Lab, use_uplab);

% -------------------------------------------------------------------------
% If dbg mode, display a figure of the outputted colormap
if dbg;
    img = repmat(cmap,[1 1 20]);
    img = permute(img,[1 3 2]);
    figure;
    imagesc(img);
    axis xy;
    title('Output colormap');
    
    plot_labcurve_rgbgamut(Lab, use_uplab);
end

end