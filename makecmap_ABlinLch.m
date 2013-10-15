function cmap = makecmap_ABlinLch(params, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<2 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end

% -------------------------------------------------------------------------
% Check input is okay
neces_fields = {'n','lch1','lch2'};
li = isfield(params,neces_fields);
if any(~li)
    error('Field %s is blank. ',neces_fields{~li});
end

% -------------------------------------------------------------------------
% Fill in non-essential fields
if ~isfield(params,'use_uplab')
    params.use_uplab = false;
end

% -------------------------------------------------------------------------
% Parse input
n         = params.n;
use_uplab = params.use_uplab;
lch1      = params.lch1;
lch2      = params.lch2;

% -------------------------------------------------------------------------
% Check lch1 and lch2 are in gamut
gamut = fetch_cielchab_gamut('srgb',[],[],use_uplab);
li = ~isingamut([lch1;lch2], gamut, 'lch');
if any(li)
    erstr = '';
    if li(1);   erstr = mat2str(lch1); end
    if all(li); erstr = [erstr ' and ']; end
    if li(2);   erstr = [erstr mat2str(lch2)]; end
    error('LCh input %s is outside srgb gamut',erstr);
end

% -------------------------------------------------------------------------
% Build the colormap

% Linear in LCh 
L = linspace(lch1(1), lch2(1), n)';
C = linspace(lch1(2), lch2(2), n)';
h = linspace(lch1(3), lch2(3), n)';

% Convert to Lab
Lab = [L C.*cosd(h) C.*sind(h)];

% Convert from Lab to srgb
cmap = hard_lab2rgb(Lab, use_uplab);

% -------------------------------------------------------------------------
% If dbg mode, display a figure of the outputted colormap
if dbg;
    % Plot the colormap
    figure;
    imagesc(permute(cmap,[1 3 2]));
    axis xy;
    title('Output colormap');
    
    % Colormap in 3d gamut
    plot_labcurve_rgbgamut(Lab, use_uplab);
end

end

% makecmap_ABlinLch(struct('n',32,'lch1',[30,50,30],'lch2',[60 5 370]),1);