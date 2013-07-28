function cmap = cie_rainbow_cmap(n, attr, verbose, func)

% -------------------------------------------------------------------------
% Default inputs
if nargin<4
    func = []; % Only used if stored cmap doesn't exist or needs replacing
end
if nargin<3 || isempty(verbose)
    verbose = 0; % Whether to output information and figures
end
if nargin<2 || isempty(attr)
    attr = 'greenmid'; % Colormap type option
end
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1); % Number of colours in the colormap
end

% -------------------------------------------------------------------------
% Parameters
n_file = 256; % Number of colours in saved cmap. 256 for smooth, 257 for greenmid.
interp_method = 'linear'; % Interpolation method used to go from saved cmap to output

% -------------------------------------------------------------------------
% Define filename for pre-parsed csv file
switch lower(attr)
    case 'smooth'
        fname = 'ciebow_cmap_smooth.tsv';
    case 'greenmid'
        fname = 'ciebow_cmap_greenmid.tsv';
    otherwise
        error('Unfamiliar colormap attribute: %s',attr);
end

dirname = fileparts(mfilename('fullpath')); % Folder containing this .m file
file = fullfile(dirname,fname);

% Try to load the pre-parsed colormap from file
if exist(file,'file')
    raw_cmap = load(file);
end

% Check this is ok
if ~exist(file,'file') || size(raw_cmap,1)<n_file
    if verbose; disp('Parsing a colormap to store for future use'); end;
    % Make a highly sampled colormap to store for future use
    raw_cmap = cie_rainbow_cmap_make(n_file, attr, func, verbose);
    % Save it as a tab deliminated file
    dlmwrite(file,raw_cmap,'delimiter','\t','precision','%.6f');
end

% Interpolate the colormap for desired n
x  = 1:size(raw_cmap,1);
xi = linspace(1,size(raw_cmap,1),n);
cmap = interp1(x,raw_cmap,xi,interp_method);

% If verbose mode, display a figure of the outputted colormap
if verbose;
    img = repmat(cmap,[1 1 20]);
    img = permute(img,[1 3 2]);
    figure;
    imagesc(img);
    axis xy;
    title('Output colormap');
end

end