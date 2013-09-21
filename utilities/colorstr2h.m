function h = colorstr2h(c,use_uplab)
    if nargin<2
        use_uplab = false;
    end
    rgb = colorstr2rgb(c);
    lab = rgb2lab(rgb,use_uplab);
    h = mod(atan2(lab(3),lab(2))/pi*180,360);
end