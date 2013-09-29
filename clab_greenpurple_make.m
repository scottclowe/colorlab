function cmap = clab_greenpurple_make(n, attr, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<3 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end
if nargin<2 || isempty(attr)
    attr = 'Lfix'; % Colormap type option
end
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1); % Number of colours in the colormap
end

% -------------------------------------------------------------------------
% Lookup parameters
params = struct('n',n);
params.use_uplab = false;

switch lower(attr)
    case 'lfix'
        %        CIE  [  L*     C*      h]
        params.lch1 = [ 66.25   94.9  328]; % purple
        params.lch2 = [ 66.25   94.9  136]; % green
    case {'lfix180','lfixhfix'}
        %        CIE  [  L*     C*      h]
        params.lch1 = [ 60.75   83    314]; % purple
        params.lch2 = [ 60.75   83    134]; % green
    case 'lvar'
        %        CIE  [  L*     C*      h]
        params.lch1 = [ 46.25  112.5  314]; % purple
        params.lch2 = [ 87.75  112.5  134]; % green
    case 'maxsep'
        %        CIE  [  L*     C*      h]
        params.lch1 = [ 30.5   130    302]; % blue
        params.lch2 = [ 90     100    122]; % green
    otherwise
        error('Unfamiliar colormap attribute: %s',attr);
end

% -------------------------------------------------------------------------
% Build the colormap
cmap = makecmap_ABlin(params, dbg);

end

% % hue found using gamut_interactive_app

% % Then Chroma and Lightness found using the following:
% k = min(g.lch(g.lch(:,3)==134,2),g.lch(g.lch(:,3)==134+180,2));
% [C,I] = max(k)
% gg = g.lch(g.lch(:,3)==134,:);
% gg(I,1)
% C