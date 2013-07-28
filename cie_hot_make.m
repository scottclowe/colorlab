function cmap = cie_hot_make(n, method, func, dbg)

if nargin<4
    dbg = 0;
end
if nargin<3
    func = [];
end
if nargin<2
    method = '1';
end

if isnumeric(method); method = num2str(method); end

% Get params
params = get_params(method);

L = linspace(params.L0,params.Lmax,n);
Lmid = (params.L0+params.Lmax)/2;
h = params.h0 + params.h_per_L * ((L-params.L0) + params.L_off - (L-Lmid).^2 * params.L_off/(params.Lmax-Lmid)^2);
c = params.maxc - abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(params.Lmax-Lmid,Lmid-params.L0))).^params.expnt) * (params.maxc-params.c0) / abs(Lmid.^params.expnt);

% Go from cylindrical Lch co-ords to cartesian Lab
a = c.*cos(h/360*(2*pi));
b = c.*sin(h/360*(2*pi));

% Turn Lab into sRGB values
Lab = [L' a' b'];
cmap = gd_lab2rgb(Lab, func);

if dbg
    % Plot the colormap
    img = repmat(cmap,[1 1 20]);
    img = permute(img,[1 3 2]);
    figure;
    imagesc(img);
    axis xy;
end

end

function params = get_params(method)

% Initial values
params.expnt    = 2;
params.h0       = 0;
params.h_per_L  = 1;
params.L_off    = 0;
params.maxc     = [];
params.L0       = 0;
params.Lmax     = 100;
params.c0       = 0;

switch method
    case '1' %25
        % Based on Set#5
        params.expnt = 2;
        params.h0 = 6.5;
        params.hmax = 106;
        params.L_off = 0;
        params.L0   = 2;
        params.Lmax = 98;
        params.maxc = 73.9; % 75.6 -> 74.25 -> 73.9
        
        params.h_per_L = (params.hmax-params.h0)/100;
        
    case '2' %9
        params.expnt = 2.2;
        params.h0 = 15;
        params.hmax = 100;
        params.L_off = 0;
        params.L0   = 3;
        params.Lmax = 97;
        params.maxc = 73; %[];
        
        params.h_per_L = (params.hmax-params.h0)/100;
        
    case '3' %11
        params.expnt = 2.3;
        params.h0 = 10;
        params.hmax = 100;
        params.L_off = 5;
        params.L0   = 3;
        params.Lmax = 97;
        params.maxc = 71; %[];
        
        params.h_per_L = (params.hmax-params.h0)/100;
        
        
    case '4' %10
        params.expnt = 2.5;
        params.h0 = 15;
        params.hmax = 100;
        params.L_off = 5;
        params.L0   = 3;
        params.Lmax = 97;
        params.maxc = 67.8; %[];
        
        params.h_per_L = (params.hmax-params.h0)/100;
    case '5' %12
        params.expnt = 2.3;
        params.h0 = 15;
        params.hmax = 100;
        params.L_off = 5;
        params.L0   = 3;
        params.Lmax = 97;
        params.maxc = 67.7; %[];
        
        params.h_per_L = (params.hmax-params.h0)/100;
        
    otherwise
        error('Unfamiliar method');

end

end