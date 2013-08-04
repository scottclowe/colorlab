function [gamut] = fetch_cielchab_gamut(space, N, point_method, use_uplab)

verbose = 0;

% -------------------------------------------------------------------------
% - INPUT HANDLING -
if nargin<1
    space = 'srgb';
end
if nargin<2
    N = [];
end
if nargin<3
    point_method = 'face';
end
if nargin<4
    use_uplab = false;
end

% -------------------------------------------------------------------------
% Swap to canonical names
switch lower(space)
    case {'srgb','rgb'}
        space = 'srgb'; % Canonical
end
switch lower(point_method)
    case {'cube','point'}
        point_method = 'cube'; % Canonical
    case {'face','faces'}
        point_method = 'face'; % Canonical
    case {'face-plus','face-cap'}
        point_method = 'face-plus'; % Canonical
end

% -------------------------------------------------------------------------
% Check to see if matfile exists
if use_uplab
    matname = [space 'gamut_uplab.mat'];
else
    matname = [space 'gamut.mat'];
end
dirname = fileparts(mfilename('fullpath')); % Folder containing this .m file
fname   = fullfile(dirname,matname);

% -------------------------------------------------------------------------
% If it doesn't, run script and save
if ~exist(fname,'file')
    if verbose; fprintf('Couldn''t load %s\n',fname); end
    [gamut] = make_cielchab_gamut(space, N, point_method, use_uplab);
    save(fname,'-struct','gamut');
    if verbose; fprintf('Made and saved gamut as %s\n',fname); end
    return;
end

% -------------------------------------------------------------------------
% If it does, load it and check settings are ok
gamut = load(fname);

% If no inputs are given, return whatever is loaded
if nargin<2; return; end

% If inputs were given, check options are okay and resolution is at least
% as good as requested
if strcmp(gamut.point_method, point_method) || ...
    (strcmp(gamut.point_method,'face-plus') && strcmp(point_method,'face'))
    if  gamut.N==N
        return;
    elseif gamut.N>N
        if verbose;
            fprintf('Returning gamut with %d sampling points per dim instead of %d\n',...
                gamut.N,N);
        end
        return;
    elseif verbose;
        fprintf('Parameters on existing matfile N=%d no good. Replacing.\n',...
            gamut.N);
    end
else
    if verbose; fprintf('Parameters on existing matfile no good. Replacing.\n'); end
end

% If options not good enough, make with given inputs and save
[gamut] = make_cielchab_gamut(space, N, point_method, use_uplab);
save(fname,'-struct','gamut');
if verbose; fprintf('Made and saved %s gamut as %s\n',space,fname); end

end