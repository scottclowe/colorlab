function rgb = makecmap_AwpBtwist(params, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<2 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end

% -------------------------------------------------------------------------
% Check input is okay
neces_fields = {'n','h1edg','h1mid','h2edg','h2mid','typ','maxc',...
    'curveLmin','curveLmax','spotLmin','spotLmax'};
li = isfield(params,neces_fields);
if any(~li)
    error('Field %s is blank. ',neces_fields{~li});
end

% -------------------------------------------------------------------------
% Fill in non-essential fields
if ~isfield(params,'use_uplab')
    params.use_uplab = false;
end
if ~isfield(params,'c0')
    params.c0 = 0;
end
% if ~isfield(params,'typ')
%     params.typ = 'sin';
% end
if ~isfield(params,'expnt')
    switch params.typ
        case 'sin'
            params.expnt = 1;
        case 'pow'
            params.expnt = 2;
    end
end

% -------------------------------------------------------------------------
% Parse input
use_uplab = params.use_uplab;
n = params.n;
h1edg = params.h1edg;
h1mid = params.h1mid;
h2edg = params.h2edg;
h2mid = params.h2mid;
curveLmin = params.curveLmin;
curveLmax = params.curveLmax;
spotLmin  = params.spotLmin;
spotLmax  = params.spotLmax;
c0    = params.c0;
expnt = params.expnt;
typ   = params.typ;
maxc  = params.maxc;

% -------------------------------------------------------------------------
% Build the colormap

% Need an odd number of colours so that grey is in the middle
neach = floor(n/2)+1;
% => if n=odd  then neach+(neach-1)=n
%    if n=even then neach+(neach-1)=n+1

L = linspace(spotLmin, spotLmax, neach);

% Chroma
curveLmid = (curveLmin+curveLmax)/2;
switch typ
    case 'sin'
        c = c0 + (1-c0) * sin(pi* (L-curveLmin)/(curveLmax-curveLmin) ).^expnt;
    case 'pow'
        c = 1 - (1-c0) * abs(((L-curveLmid)*(min(curveLmid-0,100-curveLmid)/min(curveLmax-curveLmid,curveLmid-curveLmin))).^expnt) / abs(curveLmid.^expnt);
        c = max(0,c);
    otherwise
        error('Unfamiliar type');
end

c = c*maxc;

% Hues
h1 = linspace(h1edg, h1mid, neach);
h2 = linspace(h2edg, h2mid, neach);

Lch1 = [L' c' h1'];
Lch2 = [L' c' h2']; Lch2 = flipud(Lch2);

% If min chroma is 0, then cut off the duplicated point
if c0==0
    Lch1 = Lch1(1:end-1,:);
end
lch = [Lch1; Lch2];

% Convert from lch to lab
lab = [lch(:,1) lch(:,2).*cosd(lch(:,3)) lch(:,2).*sind(lch(:,3))];
% lch = [lab(1) sqrt(lab(2)^2+lab(3)^2) mod(atan2(lab(3),lab(2))/pi*180,360)];

% Convert to rgb
rgb = soft_lab2rgb(lab, use_uplab);


% -------------------------------------------------------------------------
% If dbg mode, display a figure of the outputted colormap
if dbg;
    % Plot the colormap
    figure;
    imagesc(permute(rgb,[1 3 2]));
    axis xy;
    title('Output colormap');
    
    % Colormap in 3d gamut
    plot_labcurve_rgbgamut(lab, use_uplab);
    
    % Construction figure
    g = fetch_cielchab_gamut('srgb', [], [], use_uplab);
    li_L   = g.lchmesh.Lvec>=spotLmin & g.lchmesh.Lvec<=spotLmax;
    jointL = g.lchmesh.Lvec(li_L)';
    % jointC = min(g.lchmesh.cgrid(ismember(g.lchmesh.hvec,[h1 h2]),li_L), 1)';
    
%     figure; set(gca,'Color',[.4663 .4663 .4663]); hold on; box on;
%     plot(g.lchmesh.cgrid(g.lchmesh.hvec==h1,li_L),jointL,'b-');
%     plot(g.lchmesh.cgrid(g.lchmesh.hvec==h2,li_L),jointL,'r-');
%     plot(c,L,'ks-');
end

end