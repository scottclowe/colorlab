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

if ~use_uplab
%     edgL  = [ 15     20    25    30    35    40    45    50    55    60    65     70      75  ];
%     edgh0 = [-90.75 -85.5 -83.5 -82.5 -82.5 -82.5 -82.5 -82.5 -82.5 -83.5 -66.0 -108.75 -145  ];
%     edgh1 = [158.75 153.0 150.5 149.5 149.5 149.5 149.5 149.5 149.5 150.5 156.75 178     215  ];
%     edgc  = [ 21.9   28.3  33.9  38.8  43.1  47.3  51.5  55.8  60.0  62.8  59.2  46.4     38.4];
    
    [h0,h1,c] = arcparams_cielab(L, rgbgamut);
    
else
    % Stuff here
%     edgL  = [ 40     45    50     55     60     65    70     75   ];
%     edgh0 = [-54    -54   -53    -46.75 -55.75 -63   -77.5  -140  ];
%     edgh1 = [ 57.75  57.5  56.75  57     67.50  86.5 193.75  220  ];
%     edgc  = [ 72.5   80    87     90     75     62.9  50      40.1];
end

% -------------------------------------------------------------------------
if h1-h0>359
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

% =========================================================================
function [h0,h1,c] = arcparams_cielab(L, rgbgamut)
    % ---------------------------------------------------------------------
    if nargin<2
        rgbgamut = fetch_cielchab_gamut('srgb',[],[],use_uplab);
    end
    % ---------------------------------------------------------------------
    %  0<L<61:  Min1 at h=214, Min2 at h=93
    % 61<L<65:  Min1 at h=214, Min2 at h=14
    % At L=65:  Min1 splits at h=55
    % 75<L<100: Min1 and Min2 level (h=22, h=273)
    
    % Define our seach space for second minima
    if L<65
        h_srt =   0;
        h_end = 135;
    elseif L<71
        h_srt = -106;
        h_end =  135;
    else
        h0 = -145;
        h1 =  215;
        c  = min(rgbgamut.lchmesh.cgrid(:, rgbgamut.lchmesh.Lvec==L));
        return;
    end
    if h_srt>=0
        hli = rgbgamut.lchmesh.hvec >= h_srt & rgbgamut.lchmesh.hvec <= h_end;
    else
        hli = rgbgamut.lchmesh.hvec <= h_end | rgbgamut.lchmesh.hvec >= (h_srt+360);
    end
    % Look up chroma
    cc = rgbgamut.lchmesh.cgrid(:, rgbgamut.lchmesh.Lvec==L);
    c  = 0.99*min(cc(hli,:));
    h0 = rgbgamut.lchmesh.hvec(find(~(cc>=c),1,'last' )+1)-360;
    h1 = rgbgamut.lchmesh.hvec(find(~(cc>=c),1,'first')-1);
end

% =========================================================================
% Search code

% cc = rgbgamut.lchmesh.cgrid(:,rgbgamut.lchmesh.Lvec==30);
% figure; plot(rgbgamut.lchmesh.hvec, cc);
% % Manually pick value for c
% h0 = rgbgamut.lchmesh.hvec(find(~(cc>=c),1,'last')+1)-360
% h1 = rgbgamut.lchmesh.hvec(find(~(cc>=c),1,'first')-1)
% 
% figure;
% imagesc(rgbgamut.lchmesh.Lvec, rgbgamut.lchmesh.hvec, bsxfun(@rdivide, rgbgamut.lchmesh.cgrid, max(rgbgamut.lchmesh.cgrid)))
% colormap(clab_hot); colorbar; xlabel('L'); ylabel('h');
% figure;
% imagesc(rgbgamut.lchmesh.Lvec, rgbgamut.lchmesh.hvec, bsxfun(@rdivide, 1./rgbgamut.lchmesh.cgrid, max(1./rgbgamut.lchmesh.cgrid)))
% colormap(clab_hot); colorbar; xlabel('L'); ylabel('h');
