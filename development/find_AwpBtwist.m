% Find seperation distance for a set of different hues
% (hot-white-cold style)
clear;
% close all;

%% Parameters

use_uplab = false;
npoints = 100;     % number of points to use in test series
typ = 'sin';       % 'pow' or 'sin'
c0rel = 0;
ncurve = 100;
nLmax  = 500;

switch typ
    case 'sin'
        expnt = 1;
    case 'pow'
        expnt = 2;
end

% Best hot-twist params (in terms of looks, mostly): Take #1
%         h1edg: 316
%         h1mid: 270
%         h2edg: 37
%         h2mid: 80
%          Ledg: 22
%          Lmid: 94
%          maxc: 70.2095
%     use_uplab: 0
%             n: 32
%            c0: 0
%           typ: 'sin'
%         expnt: 1
%         Lmaxc: 43.5208

% Best hot-twist params (in terms of looks, mostly): Take #2
%         h1edg: 318
%         h1mid: 270
%         h2edg: 27
%         h2mid: 92
%          Ledg: 17
%          Lmid: 92
%          maxc: 71.2675
%     use_uplab: 0
%             n: 32
%           typ: 'sin'
%         expnt: 1
%         Lmaxc: 44.8764
%            c0: 0

% Best hot-twist params (in terms of looks, mostly): Take #3
%         h1edg: 318
%         h1mid: 270
%         h2edg: 22
%         h2mid: 92
%          Ledg: 19
%          Lmid: 88
%          maxc: 74.9355
%     use_uplab: 0
%             n: 32
%           typ: 'sin'
%         expnt: 1
%         Lmaxc: 45.6301
%            c0: 0

%% Multi-loop

g = fetch_cielchab_gamut('srgb', [], [], use_uplab);
glchmesh = g.lchmesh;

% Don't vary more than a couple of dimensions at once,
% or curse of dimensionality will bite you.
mp_names = { 'h1edg'  , 'h1mid'  , 'h2edg', 'h2mid' , 'Ledg' , 'Lmid' };
% mp_inpts = { 290:4:320, 250:4:270, 20:4:40, 80:4:100, 25, 90};
% mp_inpts = { 315:4:327, 266:4:280, 32:4:40, 76:4:84, 25, 95};
% mp_inpts = { 310:1:316, 238:2:270, 36:.25:38, 80:4:102.8, 22, 94};
% mp_inpts = { 314:1:316, 278:2:300, 36:.25:38, 71:1:82, 22, 94};
% mp_inpts = { 318, 270, 27, 92, 17, 92};
% mp_inpts = { 318, 270, 20, 92, 17, 88};
% mp_inpts = { 318, 270, 17, 95, 17, 86};
% mp_inpts = { 318, 270, 19, 95, 20, 87};
% mp_inpts = { 318, 270, 21, 95, 20, 88};
% mp_inpts = { 318, 270, 19, 95, 18, 88};
% mp_inpts = { 318, 270, 21, 92, 18, 88};
% mp_inpts = { 318, 270, 22, 92, 19, 88};
% mp_inpts = { 318, 270, 22, 92, 19, 88};
% mp_inpts = { 318, 270, 22, 92, 19, 88};
% mp_inpts = { 318, 270, 25, 90, 19, 89};
mp_inpts = { 318, 270, 23, 90, 19, 88};

mp_sizes = cellfun('prodofsize',mp_inpts);
% mp_indxs = ones(size(mp_names)); % ones(1,numel(mp_names));
% mp_indxs(1) = 0;
mp_outps = nan(mp_sizes);
mp_Lmids = nan(mp_sizes);

end_itern = prod(mp_sizes);
iteration = 0;
echo_each = 200;

if matlabpool('size')==0; matlabpool open; end;
srt = tic;

mp_sizes_cp = [1 cumprod(mp_sizes(1:end-1))];

% Simulate a stacked for-loop
% Could do this as a parfor loop
parfor iteration=1:end_itern
    
    mp_indxs = {};
    [mp_indxs{1:numel(mp_sizes)}] = ind2sub(mp_sizes,iteration);
    params = struct;
    for j=1:numel(mp_names)
        params.(mp_names{j}) = mp_inpts{j}(mp_indxs{j});
    end
    
    % ---------------------------------------------------------------------
    % Now run the thing
    L  = linspace(params.Ledg , params.Lmid , ncurve)';
    h1 = linspace(params.h1edg, params.h1mid, ncurve)';
    h2 = linspace(params.h2edg, params.h2mid, ncurve)';
    
    % Interpolate on mesh to find C for curve given by L and h in each half
    maxC1 = interp2(glchmesh.Lgrid, glchmesh.hgrid, glchmesh.cgrid, L, h1);
    maxC2 = interp2(glchmesh.Lgrid, glchmesh.hgrid, glchmesh.cgrid, L, h2);
    % We will mirror C across both halves, so find the min of the pair
    maxC  = min(maxC1,maxC2);
    
    % Generate a set of higher precision L values as candidates for
    % maximum chroma point
    % This has to be in the lower half of the curve because C is tighter in
    % the middle than the edge
    LmaxCs = linspace((params.Lmid+params.Ledg)/2, params.Ledg, nLmax);
    % This is a row vector, so we make a matrix of possibilities
    
    switch typ
        case 'sin'
            % Try a cosine curve from Ledg to Lmid with a start at Lmid and
            % a peak (pi/2) at LmaxCs(i)
            c = c0rel + (1-c0rel) * cos(pi* bsxfun(@rdivide,bsxfun(@minus,L,LmaxCs),2*abs(params.Lmid-LmaxCs)) ).^expnt;
        case 'pow'
            c = 1 - (1-c0rel) * bsxfun(@rdivide,bsxfun(@minus,L,LmaxCs),abs(params.Lmid-LmaxCs)).^expnt;
            c = max(0,c);
        otherwise
            error('Unfamiliar type');
    end
    % Divide gamut edge by generated curve C values to find maxc,
    % which is our scaling factor. Min to get best obtainable for
    % individual candidate. Max to get best candidate maxc.
    my_maxc = max(min(bsxfun(@rdivide,maxC,c)));
    
%     ndx = 1 + sum((mp_indxs(:)-1).*mp_sizes_cp(:));
%     if ndx ~= iteration; fprintf('no match\n'); end;
    mp_outps(iteration) = my_maxc;
    
    % ---------------------------------------------------------------------
    if mod(iteration,echo_each)==0
        fprintf('%d/%d in %d sec %s / %s\n',iteration,end_itern,toc(srt),mat2str(cell2mat(mp_indxs)),mat2str(mp_sizes));
    end
end

toc(srt);

%%
% Do max in each dimension and plot

clrs = clab_Lfix_disting(numel(mp_names)+1,[35 60],0);

figure; hold on;
for j=1:numel(mp_names)
    M = mmax(mp_outps,j,'x');
    plot(M(:),'-o','Color',clrs(j,:));
end
legend(mp_names,'Location','EastOutside');
title('Max');

figure; hold on;
for j=1:numel(mp_names)
    M = mmean(mp_outps,j,'x');
    plot(M(:),'-o','Color',clrs(j,:));
end
legend(mp_names,'Location','EastOutside');
title('Mean');

%%

params = struct;

[M,I] = max(mp_outps(:));
clear mp_best;
[mp_best{1:numel(mp_sizes)}] = ind2sub(mp_sizes,I);

for j=1:numel(mp_names)
    eval(sprintf('params.%s = %s;', mp_names{j}, num2str(mp_inpts{j}(mp_best{j}))));
end
params.maxc  = M;

params.use_uplab = use_uplab;
params.n     = 32;
params.typ   = typ;
params.expnt = expnt;

% Need to work out best curve shape for these parameters
L  = linspace( params.Ledg,  params.Lmid, 4*ncurve)';
h1 = linspace(params.h1edg, params.h1mid, 4*ncurve)';
h2 = linspace(params.h2edg, params.h2mid, 4*ncurve)';

maxC1 = interp2(g.lchmesh.Lgrid, g.lchmesh.hgrid, g.lchmesh.cgrid, L, h1);
maxC2 = interp2(g.lchmesh.Lgrid, g.lchmesh.hgrid, g.lchmesh.cgrid, L, h2);
maxC  = min(maxC1,maxC2);

% Increase precision on Lmaxc even further
LmaxCs = linspace((params.Lmid+params.Ledg)/2, params.Ledg, 4*nLmax);
switch typ
    case 'sin'
        % Try a cosine curve from Ledg to Lmid with a start at Lmid and
        % a peak (pi/2) at LmaxCs(i)
        C = c0rel + (1-c0rel) * cos(pi* bsxfun(@rdivide,bsxfun(@minus,L,LmaxCs),2*abs(params.Lmid-LmaxCs)) ).^params.expnt;
    case 'pow'
        C = 1 - (1-c0rel) * bsxfun(@rdivide,bsxfun(@minus,L,LmaxCs),abs(params.Lmid-LmaxCs)).^params.expnt;
        C = max(0,C);
    otherwise
        error('Unfamiliar type');
end
% Divide gamut edge by generated curve C values to find maxc,
% which is our scaling factor. Min to get best obtainable for
% individual candidate. Max to get best candidate maxc.
[my_maxc,I] = max(min(bsxfun(@rdivide,maxC,C)));
my_Lmaxc = LmaxCs(I);


params.Lmaxc = my_Lmaxc;
params.maxc  = my_maxc;
params.c0    = c0rel*params.maxc;

rgb = makecmap_AwpBtwist(params, true);
