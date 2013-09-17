function [rgb,params] = clab_Larc(n, L, use_uplab, dbg)

% -------------------------------------------------------------------------
% Default with same number of colors as in use for current colormap
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1);
end
if nargin<2 || isempty(L)
    L = 60; % Lightness
end
if nargin<3
    use_uplab = false;
end
if nargin<4
    dbg = true;
end

use_cmax = false;

% -------------------------------------------------------------------------
rgbgamut = fetch_cielchab_gamut('srgb',[],[],use_uplab);
% cc = rgbgamut.lchmesh.cgrid(:,rgbgamut.lchmesh.Lvec==L);
% figure; plot(rgbgamut.lchmesh.hvec, cc);
% % Manually pick value for c
% h0 = rgbgamut.lchmesh.hvec(find(~(cc>=c),1,'last')+1)-360
% h1 = rgbgamut.lchmesh.hvec(find(~(cc>=c),1,'first')-1)

if use_uplab
    % UPLab
    switch L
        case 40
            c = 72.5;
            h0 = -54;
            h1 =  57.75;
        case 45
            c = 80;
            h0 = -54;
            h1 =  57.5;
        case 50
%             c = 38.1;
%             h0 = -80;
%             h1 = 195;
            c = 87;
            h0 = -53;
            h1 =  56.75;
        case 55
%             c = 70;
%             h0 = -57;
%             h1 =  68.25;
            c = 90;
            h0 = -46.75;
            h1 =  57;
        case 60
%             c = 44.4;
%             h0 = -79.75;
%             h1 = 193;
            c = 75;
            h0 = -55.75;
            h1 =  67.50;
        case 65
            c = 62.9;
            h0 = -63;
            h1 =  86.5;
        case 70
            c = 50;
            h0 = -77.5;
            h1 = 193.75;
        case 75
            c = 40.1;
            h0 = -140;
            h1 =  220;
        otherwise
            h0 = -54;
            h1 =  57.5;
            if h0>=0
                hli = rgbgamut.lchmesh.hvec >= h0 & rgbgamut.lchmesh.hvec <= h1;
            else
                hli = rgbgamut.lchmesh.hvec <= h1 | rgbgamut.lchmesh.hvec >= (h0+360);
            end
            Lli = rgbgamut.lchmesh.Lvec==L;
            c = min(rgbgamut.lchmesh.cgrid(hli,Lli));
    end
else
    % CIELab
    switch L
        case 35
            c = 43.1;
            h0 = -82.5;
            h1 = 149.5;
        case 40
            c = 47.3;
            h0 = -82.5;
            h1 = 149.5;
        case 45
            c = 51.5;
            h0 = -82.5;
            h1 = 149.5;
        case 50
            c = 55.8;
            h0 = -82.5;
            h1 = 149.5;
        case 55
            c = 60.0;
            h0 = -82.5;
            h1 = 149.5;
        case 60
            c = 62.8;
            h0 = -83.5; % 276.5-360
            h1 = 150.5;
        case 65
            c = 59.2;
            h0 = -66.00;
            h1 = 156.75;
        case 70
            c = 46.4; % c = 38.4;
            h0 = -108.75;
            h1 =  178.00;
        case 75
            c = 38.4;
            h0 = -145;
            h1 =  215;
        otherwise
            h0 = -82.5;
            h1 = 149.5;
            if h0>=0
                hli = rgbgamut.lchmesh.hvec >= h0 & rgbgamut.lchmesh.hvec <= h1;
            else
                hli = rgbgamut.lchmesh.hvec <= h1 | rgbgamut.lchmesh.hvec >= (h0+360);
            end
            Lli = rgbgamut.lchmesh.Lvec==L;
            c = min(rgbgamut.lchmesh.cgrid(hli,Lli));
    end
end

% -------------------------------------------------------------------------
if h1-h0==360
    h = linspace(h0, h1, (n+1)).';
    h = h(1:end-1);
else
    h = linspace(h0, h1, n).';
end

LL = repmat(L,n,1);

if use_cmax
    P = [LL repmat(c,n,1) h];
    [TF,P2] = isingamut(P, rgbgamut, 'lch');
    c = P2(:,2);
end
a = c.*cosd(h);
b = c.*sind(h);


Lab = [LL a b];

rgb = hard_lab2rgb(Lab, use_uplab);

% -------------------------------------------------------------------------
% If dbg mode, display a figure of the outputted colormap
if dbg;
    img = repmat(rgb,[1 1 20]);
    img = permute(img,[1 3 2]);
    figure;
    imagesc(img);
    axis xy;
    title('Output colormap');
    plot_labcurve_rgbgamut(Lab, use_uplab);
end

% -------------------------------------------------------------------------
if nargout<2; return; end
params.n  = n;
params.L  = L;
params.c  = c;
params.h0 = h0;
params.h1 = h1;

end