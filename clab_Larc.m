function [varargout] = clab_Larc(n, LL, use_uplab, use_cmax, dbg)

% -------------------------------------------------------------------------
% Default with same number of colors as in use for current colormap
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1);
end
if nargin<2 || isempty(LL)
    LL = 55; % Lightness
end
if nargin<3 || isempty(use_uplab)
    use_uplab = false;
end
if nargin<4 || isempty(use_cmax)
    use_cmax = false;
end
if nargin<5
    dbg = false;
end

% -------------------------------------------------------------------------
if use_cmax && n>10
    warning('Using maximum chroma for each colour with lots of colours is a bad idea');
end

% -------------------------------------------------------------------------
rgbgamut = fetch_cielchab_gamut('srgb',[],[],use_uplab);

if ~use_uplab
%     % Use obsevered LL values only (could interpolate between these)
%     edgL  = [ 15     20    25    30    35    40    45    50    55    60    65     70      75  ];
%     edgh0 = [-90.75 -85.5 -83.5 -82.5 -82.5 -82.5 -82.5 -82.5 -82.5 -83.5 -66.0 -108.75 -145  ];
%     edgh1 = [158.75 153.0 150.5 149.5 149.5 149.5 149.5 149.5 149.5 150.5 156.75 178     215  ];
%     edgc  = [ 21.9   28.3  33.9  38.8  43.1  47.3  51.5  55.8  60.0  62.8  59.2  46.4     38.4];
    
%     % Use automated definitions based on observed search spaces
%     [h0,h1,c] = arcparams_cielab(LL, rgbgamut);
    
    % Use same h0 and h1 for all LL
    h0 = -82.5;
    h1 = 149.5;
    
else
%     % Use obsevered LL values only (could interpolate between these)
%     edgL  = [ 40     45    50     55     60     65    70     75   ];
%     edgh0 = [-54    -54   -53    -46.75 -55.75 -63   -77.5  -140  ];
%     edgh1 = [ 57.75  57.5  56.75  57     67.50  86.5 193.75  220  ];
%     edgc  = [ 72.5   80    87     90     75     62.9  50      40.1];
    
%     % Use automated definitions based on observed search spaces
%     [h0,h1,c] = arcparams_uplab(LL, rgbgamut);
    
    % Use same h0 and h1 for all LL
    h0 = -72; % include all above second-min chroma: h0 = -72;
%     h0 = -49; % start at the top of the blue ridge:  h0 = -49;
    h1 = 190;
    
end

% -------------------------------------------------------------------------
if h1-h0>359
    h = linspace(h0, h1, (n+1)).';
    h = h(1:end-1);
else
    h = linspace(h0, h1, n).';
end

if use_cmax
    cc = nan(n,length(LL));
    for i=1:length(LL)
        P = [repmat(LL(i),[n 1]) ones(n,1) h];
        [TF,P2] = isingamut(P, rgbgamut, 'lch');
        cc(:,i) = P2(:,2);
    end
else
    if h0>=0
        hli = rgbgamut.lchmesh.hvec >= h0 & rgbgamut.lchmesh.hvec <= h1;
    else
        hli = rgbgamut.lchmesh.hvec <= h1 | rgbgamut.lchmesh.hvec >= (h0+360);
    end
    % Look up chroma
    cc = min(rgbgamut.lchmesh.cgrid(hli, ismember(rgbgamut.lchmesh.Lvec,LL)));

end

aa = bsxfun(@times, cc, cosd(h));
bb = bsxfun(@times, cc, sind(h));

% -------------------------------------------------------------------------
% Now pick out these colours

LL_Lab = cell(length(LL),1);
LL_rgb = cell(length(LL),1);

for j=1:length(LL)
    L = repmat(LL(j),n,1);
    a = aa(:,j);
    b = bb(:,j);
    Lab = [L a b];
    rgb = hard_lab2rgb(Lab, use_uplab);
    LL_Lab{j} = Lab;
    LL_rgb{j} = rgb;
end

% -------------------------------------------------------------------------
% Output all in one matrix, or in individual matrices for each L
% depending on number of outputs
if nargout>1
    varargout = LL_rgb;
else
    rgb = cell2mat(LL_rgb);
    I = bsxfun(@plus,(1:n),n*(0:(length(LL)-1))');
    varargout = {rgb(I(:),:)};
end

% -------------------------------------------------------------------------
% If dbg mode, display a figure of the outputted colormap
if dbg
    rgb = cell2mat(varargout);
    figure;
    imagesc(permute(rgb,[1 3 2]));
    axis xy;
    title('Output colormap');
    
    LL_Lab2 = nan(n*length(LL),3);
    I = setdiff(1:((n+1)*length(LL)),(n+1)*(1:length(LL)));
    LL_Lab2(I,:) = cell2mat(LL_Lab);
    plot_labcurve_rgbgamut(LL_Lab2, use_uplab);
end

% -------------------------------------------------------------------------
% if nargout<2; return; end
% params.n  = n;
% params.LL  = LL;
% params.c  = c;
% params.h0 = h0;
% params.h1 = h1;

end

% =========================================================================
function [h0,h1,c] = arcparams_cielab(LL, rgbgamut)
    % ---------------------------------------------------------------------
    if nargin<2
        rgbgamut = fetch_cielchab_gamut('srgb',[],[],false);
    end
    % ---------------------------------------------------------------------
    %  0<LL<61:  Min1 at h=214, Min2 at h=93
    % 61<LL<65:  Min1 at h=214, Min2 at h=14
    % At LL=65:  Min1 splits at h=55
    % 75<LL<100: Min1 and Min2 level (h=22, h=273)
    
    % Define our seach space for second minima
    if LL<65
        h_srt =   0;
        h_end = 135;
    elseif LL<71
        h_srt = -106;
        h_end =  135;
    else
        h0 = -145;
        h1 =  215;
        c  = min(rgbgamut.lchmesh.cgrid(:, rgbgamut.lchmesh.Lvec==LL));
        return;
    end
    if h_srt>=0
        hli = rgbgamut.lchmesh.hvec >= h_srt & rgbgamut.lchmesh.hvec <= h_end;
    else
        hli = rgbgamut.lchmesh.hvec <= h_end | rgbgamut.lchmesh.hvec >= (h_srt+360);
    end
    % Look up chroma
    cc = rgbgamut.lchmesh.cgrid(:, rgbgamut.lchmesh.Lvec==LL);
    c  = 0.99*min(cc(hli,:));
    h0 = rgbgamut.lchmesh.hvec(find(~(cc>=c),1,'last' )+1)-360;
    h1 = rgbgamut.lchmesh.hvec(find(~(cc>=c),1,'first')-1);
end

% =========================================================================
function [h0,h1,c] = arcparams_uplab(LL, rgbgamut)
    % ---------------------------------------------------------------------
    if nargin<2
        rgbgamut = fetch_cielchab_gamut('srgb',[],[],true);
    end
    % ---------------------------------------------------------------------
    
    % Define our seach space for second minima
    % This was done very quickly and is not optimal
    if LL<70
        h_srt = 0;
        h_end = 166;
    else
        h0 = -145;
        h1 =  215;
        c  = min(rgbgamut.lchmesh.cgrid(:, rgbgamut.lchmesh.Lvec==LL));
        return;
    end
    
    if h_srt>=0
        hli = rgbgamut.lchmesh.hvec >= h_srt & rgbgamut.lchmesh.hvec <= h_end;
    else
        hli = rgbgamut.lchmesh.hvec <= h_end | rgbgamut.lchmesh.hvec >= (h_srt+360);
    end
    % Look up chroma
    cc = rgbgamut.lchmesh.cgrid(:, rgbgamut.lchmesh.Lvec==LL);
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
% colormap(clab_hot); colorbar; xlabel('LL'); ylabel('h');
% 
% figure;
% imagesc(rgbgamut.lchmesh.Lvec, rgbgamut.lchmesh.hvec, bsxfun(@rdivide, 1./rgbgamut.lchmesh.cgrid, max(1./rgbgamut.lchmesh.cgrid)))
% colormap(clab_hot); colorbar; xlabel('LL'); ylabel('h');
