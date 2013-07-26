% CIEDE Color difference
function dE = ciede(varargin)

% Taken from
% https://en.wikipedia.org/wiki/Color_difference#CIEDE2000
% On 2013-06-27

switch nargin
    case 2
%         if nargin==3
%             dim = varargin{3};
%         else
%             dim = 2;
%         end
        [LCh1,LCh2] = varargin{:};
        L1 = LCh1(:,1);
        C1 = LCh1(:,2);
        h1 = LCh1(:,2);
        L2 = LCh2(:,1);
        C2 = LCh2(:,2);
        h2 = LCh2(:,2);
    case 6
        [L1,C1,h1,L2,C2,h2] = varargin{:};
    otherwise
        error('Invalid number of arguments')
end


kL = 1;
kC = 1;
kH = 1;


a1 = C1.*cosd(h1);
b1 = C1.*sind(h1);
a2 = C2.*cosd(h2);
b2 = C2.*sind(h2);


dLp = L2-L1;
Lb  = (L1+L2)/2;
Cb  = (C1+C2)/2;
ap1 = a1+a1/2.*(1-sqrt(Cb.^7./(Cb.^7+25^7)));
ap2 = a2+a2/2.*(1-sqrt(Cb.^7./(Cb.*7+25^7)));
Cp1 = sqrt(ap1.^2+b1.^2);
Cp2 = sqrt(ap2.^2+b2.^2);
Cbp = (Cp1+Cp2)/2;
dCp = Cp2-Cp1;
hp1 = mod(atan2(b1,ap1)/pi*180,360);
hp2 = mod(atan2(b2,ap2)/pi*180,360);

dhp = hp2-hp1;
dhp(abs(dhp)>180 & hp2<=hp1) = dhp(abs(dhp)>180 & hp2<=hp1) + 360;
dhp(abs(dhp)>180 & hp2>hp1)  = dhp(abs(dhp)>180 & hp2>hp1)  - 360;
dhp(Cp1==0 | Cp2==0) = 0;

dHp = 2*sqrt(Cp1.*Cp2).*sind(dhp/2);
Hbp = (hp1+hp2)/2;
Hbp(abs(hp1-hp2)>180) = (hp1(abs(hp1-hp2)>180)+hp2(abs(hp1-hp2)>180)+360)/2;
Hbp(Cp1==0 | Cp2==0) = hp1(Cp1==0 | Cp2==0) + hp2(Cp1==0 | Cp2==0);

T   = 1 - 0.17*cosd(Hbp-30) + 0.24*cosd(2*Hbp) + 0.32*cosd(3*Hbp+6) - 0.20*cosd(4*Hbp-63);
SL  = 1 + 0.015*(Lb-50).^2./sqrt(20+(Lb-50).^2);
SC  = 1 + 0.045*Cbp;
SH  = 1 + 0.015*Cbp.*T;
RT  = -2*sqrt(Cbp.^7./(Cbp.^7+25^7)).*sind(60*exp(-((Hbp-275)/25).^2));

dE = sqrt((dLp/kL./SL).^2+(dCp/kC./SC).^2+(dHp/kH./SH).^2+RT.*dCp/kC./SC.*dHp/kH./SH);

end