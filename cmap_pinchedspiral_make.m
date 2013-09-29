function rgb = cmap_pinchedspiral_make(params, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<2 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end

% -------------------------------------------------------------------------
% Check input is okay
neces_fields = {'n','h1','h_per_L','Lmin','Lmax','typ','maxc'};
li = isfield(params,neces_fields);
if any(~li)
    error('Field %s is blank. ',neces_fields{~li});
end

% -------------------------------------------------------------------------
% Fill in non-essential fields
if ~isfield(params,'use_uplab')
    params.use_uplab = false;
end
if ~isfield(params,'L_off')
    params.L_off = 0;
end
if ~isfield(params,'c0')
    params.c0 = 0;
end
if ~isfield(params,'expnt')
    switch params.typ
        case 'sin'
            params.expnt = 1;
        case 'pow'
            params.expnt = 2;
    end
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
rgb = hard_lab2rgb(Lab, params.use_uplab);

% -------------------------------------------------------------------------
if dbg
    % Plot the colormap
    img = repmat(rgb,[1 1 20]);
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
    
    % Colormap in 3d gamut
    plot_labcurve_rgbgamut(Lab, params.use_uplab);
end

end