function lch = lab2lch(lab)

if size(lab,2)~=3
    error('Input must be n-by-3 matrix');
end

L = lab(:,1);
C = sqrt(lab(:,2).^2 + lab(:,3).^2);
h = mod(atan2(lab(:,3),lab(:,2))/pi*180,360);

lch = [L C h];

end