use_uplab = false;

n = 32;
h = 40;
R = 40;
Lctr = 50;

t = linspace(-pi/2,pi/2,n);
C = R*cos(t);
L = R*sin(t)+Lctr;

lch = [L' C' repmat(h,[n 1])];
lab = [lch(:,1) lch(:,2).*cosd(lch(:,3)) lch(:,2).*sind(lch(:,3))];

rgb = hard_lab2rgb(lab,use_uplab);


    % Plot the colormap
    figure;
    imagesc(permute(rgb,[1 3 2]));
    axis xy;
    title('Output colormap');
    
    % Colormap in 3d gamut
    plot_labcurve_rgbgamut(lab, use_uplab);
    