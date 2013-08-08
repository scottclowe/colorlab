function plot_labcurve_rgbgamut(Lab, use_uplab)

if nargin<2
    use_uplab = false;
end

figure;
hold on;
plot3(Lab(:,2), Lab(:,3), Lab(:,1), 'ks-');

% Get a mesh version of the gamut
rgbgamut = fetch_cielchab_gamut('srgb', [], [], use_uplab);
if ~isfield(rgbgamut,'lchmesh')
    rgbgamut.lchmesh = make_gamut_mesh(rgbgamut);
end

L = rgbgamut.lchmesh.Lgrid([1:end 1],:);
c = rgbgamut.lchmesh.cgrid([1:end 1],:);
h = rgbgamut.lchmesh.hgrid([1:end 1],:)/180*pi;
a = c.*cos(h);
b = c.*sin(h);

rgb = soft_lab2rgb([L(:) a(:) b(:)], use_uplab);

hs = surf(a,b,L,reshape(rgb,[size(L) 3]));
set(hs,'EdgeColor','none');
set(hs,'FaceAlpha',0.75);

set(gca,'Color',[0.4663 0.4663 0.4663]);
set(gca,'XLim',[-150 150],'YLim',[-150 150],'ZLim',[0 100]);
xlabel('a*')
ylabel('b*')
zlabel('L*')

view(123,10);

end