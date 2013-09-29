function varargout = clab_Lfix_disting(n, LL, greyindex, use_uplab, func, dbg)

% -------------------------------------------------------------------------
% Default with same number of colors as in use for current colormap
if nargin<1 || isempty(n)
    n = 8; %size(get(gcf,'colormap'),1);
end
if nargin<2 || isempty(LL)
    LL = 55; % Lightness
end
if nargin<3 || isempty(greyindex)
    greyindex = 1;
end
if nargin<4
    use_uplab = false;
end
if nargin<5 || isempty(func)
    func = 'min'; % 'min' | 'mean' | 'max'
end
if nargin<6
    dbg = true;
end

% -------------------------------------------------------------------------
LL = LL(:)';
if ischar(func)
    hfunc = str2func(func);
end

% -------------------------------------------------------------------------
% Find the optimally spaced hues at given L
g = fetch_cielchab_gamut('srgb',[],[],use_uplab);

cc = g.lchmesh.cgrid(:, ismember(g.lchmesh.Lvec,LL));
hh = repmat(g.lchmesh.hvec', [1 length(LL)]);
aa = cc.*cosd(hh);
bb = cc.*sind(hh);

switch func
    case {'min','max'}
        xcc = hfunc(cc,[],2);
    case 'mean'
        xcc = hfunc(cc,2);
    otherwise
        xcc = hfunc(cc);
end

if n>1
    % greyindex=0 treated like greyindex>n
    if greyindex>=0
        min_dE = xcc;
    else
        min_dE = nan(size(hh,1),1);
    end
end

if greyindex<1 || greyindex>n
    iend = n;
else
    iend = n-1;
end

pickedI = nan(iend,1);
[C, pickedI(1)] = max(xcc);

for i=2:iend
    
    % Simple 2-norm. Could swap to CIEDE2000 instead
    dE = sqrt(...
        bsxfun(@minus, aa, aa(pickedI(i-1),:)) .^2 ...
        + bsxfun(@minus, bb, bb(pickedI(i-1),:)) .^2 ...
    );
    
    switch func
        case {'min','max'}
            xdE = hfunc(dE,[],2);
        case 'mean'
            xdE = hfunc(dE,2);
        otherwise
            xdE = hfunc(dE);
    end
    min_dE = min(min_dE, xdE);
    
    [C, pickedI(i)] = max(min_dE);
    
end
% Should maximise the minimum of distances for each LL

% -------------------------------------------------------------------------
% Now pick out these colours at the specified Lightness values

LL_Lab = cell(length(LL),1);
LL_rgb = cell(length(LL),1);

for j=1:length(LL)
    c = g.lchmesh.cgrid(pickedI, g.lchmesh.Lvec==LL(j));
    h = g.lchmesh.hvec(pickedI)';
    L = repmat(LL(j),size(h));
    a = c.*cosd(h);
    b = c.*sind(h);
    Lab = [L a b];
    if greyindex>=1 && greyindex<=n
        Lab = [Lab(1:(greyindex-1),:); LL(j) 0 0; Lab(greyindex:end,:)];
    end
    rgb = hard_lab2rgb(Lab, use_uplab);
    LL_Lab{j} = Lab;
    LL_rgb{j} = rgb;
end

% -------------------------------------------------------------------------
% Output all in one matrix, or in individual matrices for each L
% depending on number of outputs
if nargout>1
    varargout = LL_rgb;
else
    rgb = cell2mat(LL_rgb);
    I = bsxfun(@plus,(1:n),n*(0:(length(LL)-1))');
    varargout = {rgb(I(:),:)};
end

% -------------------------------------------------------------------------
% If dbg mode, display a figure of the outputted colormap
if dbg;
    
    rgb = cell2mat(varargout);
    figure;
    imagesc(permute(rgb,[1 3 2]));
    axis xy;
    title('Output colormap');
    
    LL_Lab2 = nan(n*length(LL),3);
    I = setdiff(1:((n+1)*length(LL)),(n+1)*(1:length(LL)));
    LL_Lab2(I,:) = cell2mat(LL_Lab);
    plot_labcurve_rgbgamut(LL_Lab2, use_uplab);
    
    figure;
    for j=1:length(LL)
        subplot(1,length(LL),j);
        set(gca,'Color',[.4663 .4663 .4663]);
        hold on;
        plot(aa(:,j), bb(:,j), 'k-');
        plot(LL_Lab{j}(:,2), LL_Lab{j}(:,3), 'w-');
        scatter(LL_Lab{j}(:,2), LL_Lab{j}(:,3), [], LL_rgb{j}, 'filled');
        xlim([-150 150]);
        ylim([-150 150]);
        xlabel('a*');
        ylabel('b*');
        title(['L = ' num2str(LL(j))]);
    end
    
end

end