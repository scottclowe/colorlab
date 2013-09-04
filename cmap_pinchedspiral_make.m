function cmap = cmap_pinchedspiral_make(params, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<2 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end

% -------------------------------------------------------------------------
% Construct points in LCh space
L = linspace(params.Lmin, params.Lmax, params.n);
Lmid = (params.Lmin+params.Lmax)/2;
params.h_per_L = (params.h2-params.h1)/(params.Lmax-params.Lmin);
h = params.h1 + params.h_per_L * ((L-params.Lmin) + params.L_off - (L-Lmid).^2 * params.L_off/(params.Lmax-Lmid)^2);
h = mod(h,360);
switch params.typ
    case 'sin'
        c = params.c0 + (1-params.c0) * sin(pi* (L-Lmin)/(params.Lmax-Lmin) ).^params.expnt;
    case 'pow'
        c = 1 - (1-params.c0) * abs(((L-Lmid)*(min(Lmid-0,100-Lmid)/min(params.Lmax-Lmid,Lmid-params.Lmin))).^params.expnt) / abs(Lmid.^params.expnt);
        c = max(0,c);
    otherwise
        error('Unfamiliar type');
end
c = c * params.maxc;

% Go from cylindrical Lch co-ords to cartesian Lab
a = c.*cosd(h);
b = c.*sind(h);
Lab = [L' a' b'];

% Turn Lab into sRGB values
cmap = hard_lab2rgb(Lab, params.use_uplab);

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
    title(sprintf('expnt = %.3f; h1 = %.3f; h/L = %.3f; L.off = %.3f',...
        params.expnt, params.h1, params.h_per_L, params.L_off));
    xlabel('Lightness');
    legend('chroma','hue','Location','NorthWest')
    
    % Check for points out of gamut
    rgbgamut = fetch_cielchab_gamut();
    [TF,P2] = isingamut(Lab,rgbgamut,'Lab');
    fprintf('%d of %d in gamut.\n',sum(TF),length(TF));
    
    % Plot c vs maxc in gamut
    figure;
    hold on;
%     plot(c1, L, 'Color', [0 .8 0]);
    plot(c , L, 'k');
    plot(sqrt(P2(:,2).^2+P2(:,3).^2), P2(:,1), 'r');
    title(sprintf('expnt = %.3f; h1 = %.3f; h/L = %.3f; L.off = %.3f',...
        params.expnt, params.h1, params.h_per_L, params.L_off));
    xlabel('chroma');
    ylabel('Lightness');
    
    
    plot_labcurve_rgbgamut(Lab, params.use_uplab);
end

end