function cmyk = lab2cmyk(Lab, use_cielab, cmyk_space)

if nargin<2
    use_cielab = false;
end
if nargin<3
    cmyk_space = 'eu'; % It is the ISO standard...
end

% http://www.ne14design.co.uk/articles/which-cmyk-profile.htm
% http://www.color.org/chardata/drsection1.xalter
% http://www.color.org/registry/index.xalter
% http://www.eci.org/en/downloads
% http://www.jpma-net.or.jp/en/
switch lower(cmyk_space)
    case 'eu'
        icc_file = 'ISOcoated_v2_300_eci.icc'; % From ECI. ISO 12647/2-2004, v2
    case 'us'
%         icc_file = 'WebCoatedSWOP2006Grade5.icc'; % From Adobe - 2008
        icc_file = 'SWOP2006_Coated5v2.icc'; % From IDEAlliance, via ICC - 2010
    case 'japan'
        icc_file = 'JapanColor2003WebCoated.icc'; % From Adobe. Might not be the best option.
    otherwise
        error('Unfamiliar CMYK space');
end

P_cmyk = iccread(icc_file);

if use_cielab
    P_uplab = iccread('CIELab_to_UPLab.icc');
    cmyk = applycform(Lab, makecform('icc', P_uplab, P_cmyk));
else
    cmyk = applycform(Lab, makecform('CLUT', P_cmyk, 'BToA0')); % Pick one of BToA0, BToA1, BToA2
end

end