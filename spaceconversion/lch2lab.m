function lab = lch2lab(lch)

siz = size(lch);

if siz(end) ~= 3
    error('Input must have size 3 in final dimension.');
end

lch = reshape(lch, [], 3);

L = lch(:,1);
a = lch(:,2) .* cosd(lch(:,3));
b = lch(:,2) .* sind(lch(:,3));

lab = [L a b];

lab = reshape(lab, siz);

end