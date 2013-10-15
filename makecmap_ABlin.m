function cmap = makecmap_ABlin(params, dbg)

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

% Go from LCh to Lab
lab1 = [lch1(1) lch1(2).*cosd(lch1(3)) lch1(2).*sind(lch1(3))];
lab2 = [lch2(1) lch2(2).*cosd(lch2(3)) lch2(2).*sind(lch2(3))];


if mod(lch1(3)-lch2(3),360)==180
    % If two end points are directly opposite,
    % we can join them together directly
    
    % Just interpolate L, a and b between Colour 1 and Colour 2
    L = linspace(lab1(1), lab2(1), n);
    a = linspace(lab1(2), lab2(2), n);
    b = linspace(lab1(3), lab2(3), n);
    
    Lab = [L' a' b'];
    
else
    % If not, we should stay in-hue and go via a shade of grey
    
    % Find lightness of central shade of grey, using similar triangles
    L_grey = (lch1(1)*lch2(2)+lch2(1)*lch1(2))/(lch1(2)+lch2(2));
    
    % Need an odd number of colours so that grey is in the middle
    neach = floor(n/2)+1;
    % => if n=odd  then neach+(neach-1)=n
    %    if n=even then neach+(neach-1)=n+1
    
    % Colour 1 to grey
    L1 = linspace(lab1(1), L_grey, neach);
    a1 = linspace(lab1(2),      0, neach);
    b1 = linspace(lab1(3),      0, neach);
    
    % Grey to Colour 2
    L2 = linspace(L_grey, lab2(1), neach);
    a2 = linspace(     0, lab2(2), neach);
    b2 = linspace(     0, lab2(3), neach);
    
    % Merge these together
    Lab1 = [L1' a1' b1'];
    Lab2 = [L2' a2' b2'];
    
    Lab  = [Lab1(1:end-1,:); Lab2];
end

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