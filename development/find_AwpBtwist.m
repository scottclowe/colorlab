% Find seperation distance for a set of different hues
% (Blue-white-red style)
clear all;
% close all;

%% Parameters

use_uplab = false;
npoints = 100;     % number of points to use in test series
typ = 'sin';       % 'pow' or 'sin'
c0 = 0;
ncurve = 100;
nLmax  = 100;
switch typ
    case 'sin'
        expnt = 1;
    case 'pow'
        expnt = 2;
end

%% Multi-loop

g = fetch_cielchab_gamut('srgb', [], [], use_uplab);

% Don't vary more than a couple of dimensions at once,
% or curse of dimensionality will bite you.
mp_names = { 'h1edg'  , 'h1mid'  , 'h2edg', 'h2mid' , 'Ledg' , 'Lmid' };
mp_inpts = { 290:4:320, 250:4:270, 20:4:40, 80:4:100, 25, 90};

mp_sizes = cellfun('prodofsize',mp_inpts);
mp_indxs = ones(size(mp_names)); % ones(1,numel(mp_names));
mp_indxs(1) = 0;
mp_outps = nan(mp_sizes);
mp_Lmids = nan(mp_sizes);

end_itern = prod(mp_sizes);
iteration = 0;
echo_each = 200;

srt = tic;

mp_sizes_cp = [1 cumprod(mp_sizes(1:end-1))];

% Simulate a stacked for-loop
% Could do this as a parfor loop
while true
    
    iteration = iteration + 1;
    
    nextj = find(mp_indxs<mp_sizes,1);
    if isempty(nextj); break; end;
    mp_indxs(1:(nextj-1)) = 1;
    mp_indxs(nextj) = mp_indxs(nextj)+1;
    for j=1:numel(mp_names)
        eval(sprintf('%s = %s;', mp_names{j}, num2str(mp_inpts{j}(mp_indxs(j)))));
    end
    
    % ---------------------------------------------------------------------
    % Now run the thing
    L  = linspace( Ledg,  Lmid, ncurve);
    h1 = linspace(h1edg, h1mid, ncurve);
    h2 = linspace(h2edg, h2mid, ncurve);
    
    maxC1 = interp2(g.lchmesh.Lgrid, g.lchmesh.hgrid, g.lchmesh.cgrid, L, h1);
    maxC2 = interp2(g.lchmesh.Lgrid, g.lchmesh.hgrid, g.lchmesh.cgrid, L, h2);
    maxC  = min(maxC1,maxC2);
    
    LmaxCs = linspace((Lmid+Ledg)/2, Ledg, nLmax);
    my_maxc = NaN;
    switch typ
        case 'sin'
            for i=1:length(LmaxCs)
                Lmaxc = LmaxCs(i);
                c = c0 + (1-c0) * cos(pi* (L-Lmaxc)/(2*abs(Lmid-Lmaxc)) ).^expnt;
                my_maxc = max(my_maxc, min(maxC./c));
            end
        case 'pow'
            for i=1:length(LmaxCs)
                Lmaxc = LmaxCs(i);
                c = 1 - (1-c0) * ((L-Lmaxc)/(abs(Lmid-Lmaxc))).^expnt;
                c = max(0,c);
                my_maxc = max(my_maxc, min(maxC./c));
            end
        otherwise
            error('Unfamiliar type');
    end
    
%     ndx = 1 + sum((mp_indxs(:)-1).*mp_sizes_cp(:));
%     if ndx ~= iteration; fprintf('no match\n'); end;
    mp_outps(iteration) = my_maxc;
    
    % ---------------------------------------------------------------------
    if mod(iteration,echo_each)==0
        fprintf('%d/%d in %d sec %s / %s\n',iteration,end_itern,toc(srt),mat2str(mp_indxs),mat2str(mp_sizes));
    end
end

%%
% Do max in each dimension and plot


%%

params = struct;

[M,I] = max(mp_outps(:));
[mp_best{1:numel(mp_sizes)}] = ind2sub(mp_sizes,I);

for j=1:numel(mp_names)
    eval(sprintf('params.%s = %s;', mp_names{j}, num2str(mp_inpts{j}(mp_best{j}))));
end
params.maxc  = M;

params.use_uplab = use_uplab;
params.n     = 32;
params.c0    = c0;
params.typ   = typ;
params.expnt = expnt;

%
% Need to work out best curve shape for these parameters
    L  = linspace( params.Ledg,  params.Lmid, ncurve);
    h1 = linspace(params.h1edg, params.h1mid, ncurve);
    h2 = linspace(params.h2edg, params.h2mid, ncurve);
    
    maxC1 = interp2(g.lchmesh.Lgrid, g.lchmesh.hgrid, g.lchmesh.cgrid, L, h1);
    maxC2 = interp2(g.lchmesh.Lgrid, g.lchmesh.hgrid, g.lchmesh.cgrid, L, h2);
    maxC  = min(maxC1,maxC2);
    
    LmaxCs = linspace((params.Lmid+params.Ledg)/2, params.Ledg, 4*nLmax);
    my_maxc = 0;
    for i=1:length(LmaxCs)
        Lmaxc = LmaxCs(i);
        switch typ
            case 'sin'
                c = c0 + (1-c0) * cos(pi* (L-Lmaxc)/(2*abs(params.Lmid-Lmaxc)) ).^expnt;
            case 'pow'
                c = 1 - (1-c0) * ((L-Lmaxc)/(abs(params.Lmid-Lmaxc))).^expnt;
                c = max(0,c);
            otherwise
                error('Unfamiliar type');
        end
        
        this_maxc = min(maxC./c);
        if this_maxc>my_maxc
            my_maxc = this_maxc;
            my_Lmaxc = Lmaxc;
        end
    end


params.Lmaxc = my_Lmaxc;
params.maxc  = my_maxc;

rgb = makecmap_AwpBtwist(params, true);
