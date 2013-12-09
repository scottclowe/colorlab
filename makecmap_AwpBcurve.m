function rgb = makecmap_AwpBcurve(params, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<2 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end

% -------------------------------------------------------------------------
% Check input is okay
neces_fields = {'n','h1','h2','typ','maxc','Ledg','Lmid','Lmaxc'};
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
if ~isfield(params,'typ')
    params.typ = 'sin';
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
% Parse input
use_uplab = params.use_uplab;
n = params.n;
Lmaxc = params.Lmaxc;
Ledg  = params.Ledg;
Lmid  = params.Lmid;
maxc  = params.maxc;
c0rel = params.c0 / params.maxc;
h1    = params.h1;
h2    = params.h2;
typ   = params.typ;
expnt = params.expnt;

% -------------------------------------------------------------------------
% Build the colormap

% Need an odd number of colours so that grey is in the middle
neach = floor(n/2)+1;
% => if n=odd  then neach+(neach-1)=n
%    if n=even then neach+(neach-1)=n+1

L = linspace(Ledg, Lmid, neach);

switch typ
    case 'sin'
        C = c0rel + (1-c0rel) * cos(pi* (L-Lmaxc)/(2*abs(Lmid-Lmaxc)) ).^expnt;
    case 'pow'
        C = 1 - (1-c0rel) * ((L-Lmaxc)/(abs(Lmid-Lmaxc))).^expnt;
        C = max(0,C);
    otherwise
        error('Unfamiliar type');
end

% Check for points out of gamut or restrict to gamut
% ... can't do this because we would need to generate a C point for all the
% ... L values in the gamut, which we haven't done

C = C*maxc;

Lch1 = [L' C' repmat(h1,neach,1)];
Lch2 = [flipud(L') flipud(C') repmat(h2,neach,1)];

% If min chroma is 0, then cut off the duplicated point
if c0rel==0
    Lch1 = Lch1(1:end-1,:);
end
lch = [Lch1; Lch2];

% Convert from lch to lab
lab = [lch(:,1) lch(:,2).*cosd(lch(:,3)) lch(:,2).*sind(lch(:,3))];
% lch = [lab(1) sqrt(lab(2)^2+lab(3)^2) mod(atan2(lab(3),lab(2))/pi*180,360)];

% Convert to rgb
rgb = hard_lab2rgb(lab, use_uplab);


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
    li_L   = g.lchmesh.Lvec>=min(Ledg,Lmid) & g.lchmesh.Lvec<=max(Ledg,Lmid);
    jointL = g.lchmesh.Lvec(li_L)';
    % jointC = min(g.lchmesh.cgrid(ismember(g.lchmesh.hvec,[h1 h2]),li_L), 1)';
    
    figure; set(gca,'Color',[.4663 .4663 .4663]); hold on; box on;
    plot(g.lchmesh.cgrid(g.lchmesh.hvec==params.h1,li_L),jointL,'b-');
    plot(g.lchmesh.cgrid(g.lchmesh.hvec==params.h2,li_L),jointL,'r-');
    plot(C,L,'ks-');
end

end