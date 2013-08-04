function cmap = cie_hot_make(n, attr, spacefun, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<4 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end
if nargin<3
    spacefun = []; % spacefuntion to map from cielab to srgb
end
if nargin<2 || isempty(attr)
    attr = '1'; % Colormap type option
end
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1); % Number of colours in the colormap
end

if isnumeric(attr); attr = num2str(attr); end

% -------------------------------------------------------------------------
% Get params
params = get_params(attr);

% -------------------------------------------------------------------------
% Construct points in LCh space
L = linspace(params.L0,params.Lmax,n);
Lmid = (params.L0+params.Lmax)/2;
h = params.h0 + params.h_per_L * ((L-params.L0) + params.L_off - (L-Lmid).^2 * params.L_off/(params.Lmax-Lmid)^2);
c = params.maxc - abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(params.Lmax-Lmid,Lmid-params.L0))).^params.expnt) * (params.maxc-params.c0) / abs(Lmid.^params.expnt);

% Go from cylindrical Lch co-ords to cartesian Lab
a = c.*cos(h/360*(2*pi));
b = c.*sin(h/360*(2*pi));

% Turn Lab into sRGB values
Lab = [L' a' b'];
cmap = gd_lab2rgb(Lab, params.use_uplab, spacefun);

% -------------------------------------------------------------------------
if dbg
    % Plot the colormap
    img = repmat(cmap,[1 1 20]);
    img = permute(img,[1 3 2]);
    figure;
    imagesc(img);
    axis xy;
    
    % Plot changes in c and h vs L
    figure;
    hold on;
    plot(L,c,'k');
    plot(L,h,'r');
    title(sprintf('expnt = %.3f; h0 = %.3f; h/L = %.3f; L.off = %.3f',...
        params.expnt, params.h0, params.h_per_L, params.L_off));
    xlabel('Lightness');
    legend('chroma','hue','Location','NorthWest')
    
    % Check for points out of gamut
    rgbgamut = fetch_cielchab_gamut();
    [TF,P2] = isingamut(Lab,rgbgamut,'Lab');
    fprintf('attr #%s has %d of %d in gamut.\n',attr,sum(TF),length(TF));
    
    % Plot c vs maxc in gamut
    figure;
    hold on;
%     plot(c1, L, 'Color', [0 .8 0]);
    plot(c , L, 'k');
    plot(sqrt(P2(:,2).^2+P2(:,3).^2), P2(:,1), 'r');
    title(sprintf('expnt = %.3f; h0 = %.3f; h/L = %.3f; L.off = %.3f',...
        params.expnt, params.h0, params.h_per_L, params.L_off));
    xlabel('chroma');
    ylabel('Lightness');
    
    
    plot_labcurve_rgbgamut(Lab)
end

end

function params = get_params(attr)

% Initial values
params.use_uplab = false;
params.expnt    = 2;
params.h0       = 0;
params.h_per_L  = 1;
params.L_off    = 0;
params.maxc     = [];
params.L0       = 0;
params.Lmax     = 100;
params.c0       = 0;

switch attr
    case '1' %25
        % Based on Set#5
        params.expnt  = 2;
        params.h0     = 6.5;
        params.hmax   = 106;
        params.L_off  = 0;
        params.L0     = 5; %2;
        params.Lmax   = 95; %98;
        params.maxc   = 73.9; % 75.6 -> 74.25 -> 73.9
        
        params.h_per_L = (params.hmax-params.h0)/100;
        
    case '2' %9
        params.expnt  = 2.2;
        params.h0     = 15;
        params.hmax   = 100;
        params.L_off  = 0;
        params.L0     = 5; %3;
        params.Lmax   = 95; %97;
        params.maxc   = 73; %[];
        
        params.h_per_L = (params.hmax-params.h0)/100;
        
    case '3' %11
        params.expnt  = 2.3;
        params.h0     = 10;
        params.hmax   = 100;
        params.L_off  = 5;
        params.L0     = 5; %3;
        params.Lmax   = 95; %97;
        params.maxc   = 71; %[];
        
        params.h_per_L = (params.hmax-params.h0)/100;
        
        
    case '4' %10
        params.expnt  = 2.5;
        params.h0     = 15;
        params.hmax   = 100;
        params.L_off  = 5;
        params.L0     = 5; %3;
        params.Lmax   = 95; %97;
        params.maxc   = 67.8; %[];
        
        params.h_per_L = (params.hmax-params.h0)/100;
        
        
    case '5' %12
        params.expnt  = 2.3;
        params.h0     = 15;
        params.hmax   = 100;
        params.L_off  = 5;
        params.L0     = 5; %3;
        params.Lmax   = 95; %97;
        params.maxc   = 67.7; %[];
        
        params.h_per_L = (params.hmax-params.h0)/100;
        
        
    case 'cool:1' %104
        params.expnt  = 2;
        params.h0     = 340;
        params.hmax   = 172;
        params.L_off  = -15;
        params.L0     = 3;
        params.Lmax   = 85;
        params.maxc   = 70; %[];
        
        params.h_per_L = (params.hmax-params.h0)/100;
        
    otherwise
        error('Unfamiliar attribute');

end

end