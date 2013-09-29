function rgb = cmap_AwpBlin_make(params, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<2 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end

% -------------------------------------------------------------------------
% Check input is okay
neces_fields = {'n','h1','h2','Ledg','Lmid','maxc'};
li = isfield(params,neces_fields);
if any(~li)
    error('Field %s is blank. ',neces_fields{~li});
end

% -------------------------------------------------------------------------
% Fill in non-essential fields
if ~isfield(params,'use_uplab')
    params.use_uplab = false;
end
if ~isfield(params,'via_black')
    params.via_black = false;
end

% -------------------------------------------------------------------------
% Parse input
use_uplab = params.use_uplab;
n  = params.n;
h1 = params.h1;
h2 = params.h2;
Ledg = params.Ledg;
Lmid = params.Lmid;
% via_black = params.via_black;
maxc  = params.maxc;

% -------------------------------------------------------------------------
% Build the colormap

% Need an odd number of colours so that grey is in the middle
neach = floor(n/2)+1;
% => if n=odd  then neach+(neach-1)=n
%    if n=even then neach+(neach-1)=n+1

L = linspace(Ledg, Lmid, neach);
c = linspace(maxc, 0, neach);
% if via_black; L = fliplr(L); end;

Lch1 = [L' c' repmat(h1,neach,1)];
Lch2 = [flipud(L') flipud(c') repmat(h2,neach,1)];

% Cut off duplicated white-point
Lch1 = Lch1(1:end-1,:);

% If n is even, cut out the white-point completely to make an even output
if mod(n,2)==0
    Lch2 = Lch2(2:end,:);
end

% Convert from lch to lab
lch = [Lch1; Lch2];
lab = [lch(:,1) lch(:,2).*cosd(lch(:,3)) lch(:,2).*sind(lch(:,3))];
% lch = [lab(1) sqrt(lab(2)^2+lab(3)^2) mod(atan2(lab(3),lab(2))/pi*180,360)];

% Convert from Lab to srgb
rgb = hard_lab2rgb(lab, use_uplab);

% -------------------------------------------------------------------------
% If dbg mode, display a figure of the outputted colormap
if dbg
    
    % Plot the colormap
    figure;
    imagesc(permute(rgb,[1 3 2]));
    axis xy;
    title('Output colormap');
    
    % Colormap in 3d gamut
    plot_labcurve_rgbgamut(lab, use_uplab);
    
    % Construction figures
    rgbgamut = fetch_cielchab_gamut('srgb', [], [], use_uplab);
    
    figure; set(gca,'Color',[.4663 .4663 .4663]); hold on; box on;
    gc1 = rgbgamut.lchmesh.cgrid(rgbgamut.lchmesh.hvec==h1, :);
    gc2 = rgbgamut.lchmesh.cgrid(rgbgamut.lchmesh.hvec==h2, :);
    plot(gc1, rgbgamut.lchmesh.Lvec, 'b-');
    plot(gc2, rgbgamut.lchmesh.Lvec, 'r-');
    plot([0 maxc],[Lmid Ledg],'k-');
    
    % Plot lines ontop of image
    figure; set(gca,'Color',[.4663 .4663 .4663]); hold on; box on;
    
    intv = 0.5;
    LL = 0:intv:100;
    cc = 0:intv:150;
    c1 = [ fliplr(cc) cc(2:end)];
    c2 = [-fliplr(cc) cc(2:end)];
    
    ccc = meshgrid(c1,LL);
    LLL = meshgrid(LL,c1)';
    hhh = [...
        repmat(h1, [length(LL) length(cc)-1]) ...
        zeros(length(LL), 1) ...
        repmat(h2 ,[length(LL) length(cc)-1]) ...
        ];
    
    aaa = ccc.*cosd(hhh);
    bbb = ccc.*sind(hhh);
    imLab = cat(3,LLL,aaa,bbb);
    
    imrgb = hard_lab2rgb(imLab,use_uplab);
    li = any(imrgb<0, 3) | any(imrgb>1, 3);
    imrgb(repmat(li,[1 1 3])) = .4663;
    
    image(c2,LL,imrgb);
    set(gca,'YDir','normal');
    xlim([-150 150]);
    ylim([0 100]);
    set(gca,'XTick',get(gca,'XTick'));
    set(gca,'XTickLabel',abs(get(gca,'XTick')));
    
    plot([-lch(1:neach,2); lch((end-neach):end,2)], ...
         [ lch(1:neach,1); lch((end-neach):end,1)], 'ks-');

end

end