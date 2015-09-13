function varargout = clab_Lfix_disting(n, LL, greyindex, Cmatch, use_uplab, func, dbg)

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
if nargin<4 || isempty(Cmatch)
    Cmatch = 0; % 0 | 1 | 2
end
if nargin<5 || isempty(use_uplab)
    use_uplab = false;
end
if nargin<6 || isempty(func)
    func = 'min'; % 'min' | 'mean' | 'max'
end
if nargin<7 || isempty(dbg)
    dbg = false;
end

% -------------------------------------------------------------------------
LL = LL(:)';
if ischar(func)
    hfunc = str2func(func);
elseif isa(func,'function_handle')
    hfunc = func;
    func = func2str(hfunc);
else
    error('Comparitor is not a function');
end

% -------------------------------------------------------------------------
% Find the optimally spaced hues at given L

% Fetch the gamut for the colour space we are using
g = fetch_cielchab_gamut('srgb',[],[],use_uplab);

% Round LL to the nearest spacing used in g
LL = round(LL/g.Lintv)*g.Lintv;

% Reduce this down to just the Lightness values in LL
cc = g.lchmesh.cgrid(:, ismember(g.lchmesh.Lvec,LL));
hh = repmat(g.lchmesh.hvec', [1 length(LL)]);
aa = cc.*cosd(hh);
bb = cc.*sind(hh);

% Apply our cross-L comparitor to optimise chroma across LL
switch func
    case {'min','max'}
        xcc = hfunc(cc,[],2);
    case 'mean'
        xcc = hfunc(cc,2);
    otherwise
        xcc = hfunc(cc);
end

% Initialise a distance from last colour measure if we are using grey
if n>1
    % greyindex=0 treated exactly like greyindex>n
    if greyindex>=0 || Cmatch
        min_dE = xcc;
    else
        min_dE = nan(size(hh,1),1);
    end
end

% If we are including grey, we need to find one fewer colour
if greyindex<1 || greyindex>n
    iend = n;
else
    iend = n-1;
end

% Pick our first colour so it has maximum Chroma
pickedI = nan(iend,1);
[C, pickedI(1)] = max(xcc);

% Loop until we have all the colours we need to find
for i=2:iend
    
    % Simple 2-norm. Could swap to CIEDE2000 instead
    % We do the 2-norm from the last colour
    % We already did all the distance from our other colours in the
    % previous runs through the loop
    dE = sqrt(...
        bsxfun(@minus, aa, aa(pickedI(i-1),:)) .^2 ...
        + bsxfun(@minus, bb, bb(pickedI(i-1),:)) .^2 ...
    );
    
    % Collapse down the different values for LL
    switch func
        case {'min','max'}
            xdE = hfunc(dE,[],2);
        case 'mean'
            xdE = hfunc(dE,2);
        otherwise
            xdE = hfunc(dE);
    end
    % Minimum of this and previous run
    min_dE = min(min_dE, xdE);
    
    % Maximise the metric to extract the new colour's index
    [C, pickedI(i)] = max(min_dE);
    
end
% Should maximise the minimum of distances for each LL

% -------------------------------------------------------------------------
% Now pick out these colours at the specified Lightness values

% Store a matrix for each L value
LL_Lab = cell(length(LL),1);
LL_rgb = cell(length(LL),1);

for j=1:length(LL)
    % Pick out all the chromas for this L
    c = g.lchmesh.cgrid(pickedI, g.lchmesh.Lvec==LL(j));
    
    % Match Chroma if required
    if Cmatch==1;
        c = min(c);
    elseif Cmatch==2;
        allc = g.lchmesh.cgrid(pickedI, ismember(g.lchmesh.Lvec,LL));
        c = min(allc(:));
    elseif Cmatch~=0
        error('Unfamiliar setting for chroma matching');
    end
    
    h = g.lchmesh.hvec(pickedI)';
    L = repmat(LL(j),size(h));
    a = c.*cosd(h);
    b = c.*sind(h);
    Lab = [L a b];
    % Insert the right shade of grey
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
    % If we want more than one output, give each matrix separately
    varargout = LL_rgb;
else
    % If we only want one (or did not ask for any) output, give all the
    % colours in one matrix, with alternating L
    rgb = cell2mat(LL_rgb); % Collapse cell to a matrix
    % Now get indices to reorder the matrix rows so they alternate
    I = bsxfun(@plus, (1:n), n*(0:(length(LL)-1))');
    varargout = {rgb(I(:),:)};
end

% -------------------------------------------------------------------------
% If debug mode, display a figure of the outputted colormap
if dbg;
    
    rgb = cell2mat(varargout);
    figure;
    imagesc(permute(rgb,[1 3 2]));
    axis xy;
    title('Output colormap');
    
    % Read the Lab values but leave a NaN between each LL
    LL_Lab2 = nan((n+1)*length(LL),3);
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