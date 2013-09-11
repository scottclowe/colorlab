function cmap = clab_blackwhite_make(n, dbg)

if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1);
end
if nargin<2
    dbg = 0;
end

use_uplab = false;

L = linspace(0,100,n)';
a = zeros(n,1);
b = zeros(n,1);

Lab = [L a b];

cmap = soft_lab2rgb(Lab, use_uplab); % Needs to be soft as white is outside by 3.8011e-05

% -------------------------------------------------------------------------
% If dbg mode, display a figure of the outputted colormap
if dbg;
    img = repmat(cmap,[1 1 20]);
    img = permute(img,[1 3 2]);
    figure;
    imagesc(img);
    axis xy;
    title('Output colormap');
    
    plot_labcurve_rgbgamut(Lab);
end

end