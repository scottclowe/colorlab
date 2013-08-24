% Outdated. see cie_hot_test2

use_uplab = false;

method = 2;

g = fetch_cielchab_gamut();

switch method
    case 1
        % Method 1:
        %          [  L  c  h  ]
        % Red max chroma
        maxc_red = [54.5 106 40];
        % Yellow max chroma
        maxc_yel = [95.75 94 106];
        % Picked by using cietools to browse with chroma setting

        h_per_L = (maxc_yel(3)-maxc_red(3)) / (maxc_yel(1)-maxc_red(1));
        h0 = maxc_red(3)-h_per_L*maxc_red(1);
        
        expnt = 2;
        % % Outcome Programmatically:
        % h0 = -47.2000;
        % h_per_L = 1.6000;
        % maxc = 75.85;

    case 2
        % Method 2: Pick values and play with them
        h0 = 0;
        h_per_L = 1.15;
        expnt = 2;
        
        % % Outcome manually
        % expnt = 2;
        % h0 = 0;
        % h_per_L = 1.15;
        % maxc = 85.66
end

L = 0:100;
h = h0 + L*h_per_L;
h = round(h);

c1 = nan(size(L));
for i=1:length(L)
    tmp = g.lch(g.lch(:,1)==L(i),:);
    [ignore,I] = min(abs(tmp(:,3)-h(i)));
    c1(i) = tmp(I,2);
end

maxc = min(c1(45:55));


Lmid = 50;
Lmax = 50;
c = maxc - abs((L-Lmid).^expnt+10*L) * maxc / abs(Lmid.^expnt+10*Lmid);

(L-Lmid

A2 * (100-lmid +A1) ^2 = maxc
A2 * (    lmid +A1) ^2 = maxc

(100-lmid+A1) ^2 = maxc/A2
(    lmid+A1) ^2 = maxc/A2

100-lmid+A1 = sqrt(maxc/A2)
    lmid+A1 = sqrt(maxc/A2)

100+2*A1 = 2*sqrt(maxc/A2)

a = c.*cos(h/360*(2*pi));
b = c.*sin(h/360*(2*pi));

Lab = [L' a' b'];
rgb = gd_lab2rgb(Lab, use_uplab);

clrs = repmat(rgb,[1 1 20]);
clrs = permute(clrs,[1 3 2]);

figure;
hold on;
plot(c1, L, 'g');
plot(c , L, 'k');
title(sprintf('h0 = %.3f; h/L = %.3f; expnt = %.3f',h0,h_per_L,expnt));

figure;
imagesc(clrs);
axis xy;
title(sprintf('h0 = %.3f; h/L = %.3f; expnt = %.3f',h0,h_per_L,expnt));
