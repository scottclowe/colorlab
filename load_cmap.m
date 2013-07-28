function cmap = load_cmap(n, makefun, attr, spacefun, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<4 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end
if nargin<3
    spacefun = []; % function to map from cielab to srgb
end
if nargin<2 || isempty(attr)
    attr = []; % Colormap type option
end

% -------------------------------------------------------------------------
% Parameters
n_file = 256; % Number of colours in saved cmap. (will be odd for divergent)
interp_method = 'linear'; % Interpolation method used to go from saved cmap to output
cmap_foldername = 'cmaps';

% -------------------------------------------------------------------------

dirname = fileparts(mfilename('fullpath')); % Folder containing this .m file
dirname = fullfile(dirname, cmap_foldername); % Folder designated for stored cmaps

funname = func2str(makefun);
if strncmp(funname(end-4:end),'_make',5); funname = funname(1:end-5); end
if isempty(attr)
    fname = [funname '.tsv'];
else
    fname = [funname '_' attr '.tsv'];
end

file = fullfile(dirname,fname);

% Try to load the pre-parsed colormap from file
if exist(file,'file')
    raw_cmap = load(file);
end

% Check this is ok
if ~exist(file,'file') || size(raw_cmap,1)<n_file
    if dbg; disp('Parsing a colormap to store for future use'); end;
    
    % Make directory if non-existent
    if ~exist(dirname,'dir'); mkdir(dirname); end
    
    % Make a highly sampled colormap to store for future use
    raw_cmap = makefun(n_file, attr, spacefun, dbg);
    
    % Save it as a tab deliminated file
    dlmwrite(file,raw_cmap,'delimiter','\t','precision','%.6f');
end

% Interpolate the colormap for desired n
x  = 1:size(raw_cmap,1);
xi = linspace(1,size(raw_cmap,1),n);
cmap = interp1(x,raw_cmap,xi,interp_method);

% If dbg mode, display a figure of the outputted colormap
if dbg;
    img = repmat(cmap,[1 1 20]);
    img = permute(img,[1 3 2]);
    figure;
    imagesc(img);
    axis xy;
    title('Output colormap');
end

end