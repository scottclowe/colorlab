function rgb = clab_bluewhitered_make(n, attr, dbg)

% -------------------------------------------------------------------------
% Default inputs
if nargin<3 || isempty(dbg)
    dbg = 0; % Whether to output information and figures
end
if nargin<2 || isempty(attr)
    attr = '';
end
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1); % Number of colours in the colormap
end

% -------------------------------------------------------------------------
% Parameters
params = struct('n',n);

switch attr
    case ''
        params.h1   = 296;
        params.h2   =  40;
        params.Ledg =  46.375;
        params.Lmid =  97.411;
        params.maxc =  93.9075;
        params.use_uplab = false;
        
    case 'lfix'
        params.h1   = 296;
        params.h2   =  40;
        params.Ledg =  46.3125;
        params.Lmid =  46.3125;
        params.maxc =  94.0625;
        params.use_uplab = false;
        
%     case 'lfixdark'
% % CIELCH      [  L*      c    h]
% %   lchblue   = [ 38.50  106.64  294];
% %   lchred    = [ 38.50   82.89   41];
% %   wp        = [ 38.50,   0,      0];
%   lchblue   = [ 38.50  105  294];
%   lchred    = [ 38.50   82   41];
%   wp        = [ 38.50,   0,   0];
        
    case 'uplab'
        params.h1   = 309;
        params.h2   =  54.75;
        params.Ledg =  50.25;
        params.Lmid =  92.777;
        params.maxc =  95.6875;
        params.use_uplab = true;
        
    case 'uplab:lfix'
        params.h1   = 309;
        params.h2   =  54.75;
        params.Ledg =  50.625;
        params.Lmid =  50.625;
        params.maxc =  96.3438;
        params.use_uplab = true;
        
    otherwise
        error('Unfamiliar colormap attribute: %s',attr);
end

% -------------------------------------------------------------------------
% Build the colormap
rgb = makecmap_AwpBlin(params, dbg);

end