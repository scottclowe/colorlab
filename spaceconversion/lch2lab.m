function lab = lch2lab(lch)

if size(lch,2)~=3
    error('Input must be n-by-3 matrix');
end

L = lch(:,1);
a = lch(:,2) .* cosd(lch(:,3));
b = lch(:,2) .* sind(lch(:,3));

lab = [L a b];

end