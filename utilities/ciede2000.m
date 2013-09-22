% CIEDE Color difference
function dE = ciede2000(Lab1,Lab2)

% Taken from
% https://en.wikipedia.org/wiki/Color_difference#CIEDE2000
% On 2013-06-27

% -------------------------------------------------------------------------
% Parameters
kL = 1;
kC = 1;
kH = 1;

% -------------------------------------------------------------------------
% Parse inputs
L1 = Lab1(:,1);
L2 = Lab2(:,1);
a1 = Lab1(:,2);
a2 = Lab2(:,2);
b1 = Lab1(:,3);
b2 = Lab2(:,3);
C1 = sqrt(a1.^2+b1.^2);
C2 = sqrt(a2.^2+b2.^2);
h1 = mod(atan2(b1,a1)/pi*180, 360);
h2 = mod(atan2(b2,a2)/pi*180, 360);

% -------------------------------------------------------------------------

% Simple difference and mean
dLp = L2-L1;
Lb  = (L1+L2)/2;
Cb  = (C1+C2)/2;

% Transform a* values
ap1 = a1+a1/2.*(1-sqrt(Cb.^7./(Cb.^7+25^7)));
ap2 = a2+a2/2.*(1-sqrt(Cb.^7./(Cb.^7+25^7)));

% Recompute C and h with new a*
Cp1 = sqrt(ap1.^2+b1.^2);
Cp2 = sqrt(ap2.^2+b2.^2);
Cbp = (Cp1+Cp2)/2;
dCp = Cp2-Cp1;
hp1 = mod(atan2(b1,ap1)/pi*180,360);
hp2 = mod(atan2(b2,ap2)/pi*180,360);

% Difference in h
dhp = hp2-hp1;
dhp(abs(dhp)>180 & hp2<=hp1) = dhp(abs(dhp)>180 & hp2<=hp1) + 360;
dhp(abs(dhp)>180 & hp2>hp1)  = dhp(abs(dhp)>180 & hp2>hp1)  - 360;
dhp(Cp1==0 | Cp2==0) = 0;

dHp = 2*sqrt(Cp1.*Cp2).*sind(dhp/2);
Hbp = (hp1+hp2)/2;
Hbp(abs(hp1-hp2)>180) = (hp1(abs(hp1-hp2)>180)+hp2(abs(hp1-hp2)>180)+360)/2;
Hbp(Cp1==0 | Cp2==0) = hp1(Cp1==0 | Cp2==0) + hp2(Cp1==0 | Cp2==0);

% Compensation weights and hue rotation term
T   = 1 - 0.17*cosd(Hbp-30) + 0.24*cosd(2*Hbp) + 0.32*cosd(3*Hbp+6) - 0.20*cosd(4*Hbp-63);
SL  = 1 + 0.015*(Lb-50).^2./sqrt(20+(Lb-50).^2);
SC  = 1 + 0.045*Cbp;
SH  = 1 + 0.015*Cbp.*T;
RT  = -2*sqrt(Cbp.^7./(Cbp.^7+25^7)).*sind(60*exp(-((Hbp-275)/25).^2));

% Weighted output
dE = sqrt((dLp/kL./SL).^2+(dCp/kC./SC).^2+(dHp/kH./SH).^2+RT.*dCp/kC./SC.*dHp/kH./SH);

end