function varargout = clab_Lfix_disting(n, LL, pickL, use_uplab, dbg)

% -------------------------------------------------------------------------
% Default with same number of colors as in use for current colormap
if nargin<1 || isempty(n)
    n = 8; %size(get(gcf,'colormap'),1);
end
if nargin<2 || isempty(LL)
    LL = 50; % Lightness
end
if nargin<3 || isempty(pickL)
    pickL = LL(1); % Lightness
end
if nargin<4
    use_uplab = false;
end
if nargin<5
    dbg = true;
end

% -------------------------------------------------------------------------
% Find the optimally spaced hues at given L
g = fetch_cielchab_gamut('srgb',[],[],use_uplab);

c = g.lchmesh.cgrid(:,g.lchmesh.Lvec==pickL);
h = g.lchmesh.hvec';
L = repmat(pickL,size(h));
a = c.*cosd(h);
b = c.*sind(h);
pickLab = [L a b];

pickedI = nan(n,1);
[c1,I] = max(c);
pickedI(1) = I;
lastpickedLab = pickLab(I,:);

if n>1
    minpickLab_dE2 = nan(size(h));
end

for i=2:n
    % Don't bother to sqrt as it is a monotonic function
    lastpickLab_dE2 = sum(bsxfun(@minus, pickLab, lastpickedLab).^2, 2);
    minpickLab_dE2 = min(minpickLab_dE2, lastpickLab_dE2);
    
    [C,I] = max(minpickLab_dE2);
    
    pickedI(i) = I;
    lastpickedLab = pickLab(I,:);
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
    img = repmat(rgb,[1 1 20]);
    img = permute(img,[1 3 2]);
    figure;
    imagesc(img);
    axis xy;
    title('Output colormap');
    plot_labcurve_rgbgamut(cell2mat(LL_Lab), use_uplab);
end

end