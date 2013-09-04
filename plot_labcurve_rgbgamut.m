function plot_labcurve_rgbgamut(Lab, use_uplab)

% -------------------------------------------------------------------------
% Default inputs
if nargin<2
    use_uplab = false;
end
if nargin<1
    Lab = [];
end

% -------------------------------------------------------------------------
% Parameters
target_Lintv = 1;
target_hintv = 1;

% -------------------------------------------------------------------------
% Get a mesh version of the gamut
rgbgamut = fetch_cielchab_gamut('srgb', [], [], use_uplab);
if ~isfield(rgbgamut,'lchmesh')
    rgbgamut.lchmesh = make_gamut_mesh(rgbgamut);
end

indIntvL = max(1,floor(target_Lintv/rgbgamut.Lintv));
indIntvh = max(1,floor(target_hintv/rgbgamut.hintv));

figure;
hold on;
if ~isempty(Lab);
    plot3(Lab(:,2), Lab(:,3), Lab(:,1), 'ks');
    plot3(Lab(:,2), Lab(:,3), Lab(:,1), 'w-');
end

L = rgbgamut.lchmesh.Lgrid([1:indIntvh:end 1],1:indIntvL:end);
c = rgbgamut.lchmesh.cgrid([1:indIntvh:end 1],1:indIntvL:end);
h = rgbgamut.lchmesh.hgrid([1:indIntvh:end 1],1:indIntvL:end);
a = c.*cosd(h);
b = c.*sind(h);

rgb = soft_lab2rgb([L(:) a(:) b(:)], use_uplab);

hs = surf(a,b,L,reshape(rgb,[size(L) 3]));
set(hs,'EdgeColor','none');
if ~isempty(Lab);
    set(hs,'FaceAlpha',0.75);
end

set(gca,'Color',[0.4663 0.4663 0.4663]);
set(gca,'XLim',[-150 150],'YLim',[-150 150],'ZLim',[0 100]);
xlabel('a*')
ylabel('b*')
zlabel('L*')

view(123,10);

end