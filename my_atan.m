function t = my_atan(x,y)

t = atan(y./x); % In radians
t(isnan(t)) = 0;
s = sign(x);
s(s==1) = 0;
t = t-pi*s;
t = mod(t,2*pi);

end