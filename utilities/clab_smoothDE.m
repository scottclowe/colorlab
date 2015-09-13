function val = clab_smoothDE(val,n,space,dEmode)

if nargin<3 || isempty(space)
    space = 'lab';
end
if nargin<4 || isempty(dEmode)
    dEmode = 'cie2000';
end

switch lower(space)
    case 'lab'
        lab = val;
    case 'lch'
        lab = lch2lab(val);
    otherwise
        lab = clab_colorspace([space '->lab'],val);
end

switch lower(dEmode)
    case {'cie1976','cie76','1976','76'}
        lab_dif = diff(lab,1,1);
        dE = sqrt(sum(lab_dif.^2,2));
    case {'cie1994','cie94','1994','94'}
        error('The CIE94 colour difference is unsupported. Please use CIE2000 which supercedes it');
    case {'cie2000','cie00','2000','00'}
        dE = ciede2000(lab(1:end-1,:),lab(2:end,:));
    case {'cmcl:c','cmc','l:c'}
        error('CMC L:C colour difference is unsupported. Suggest you use CIE2000 instead.');
    otherwise
        error('Unfamiliar color distance method: %s',dEmode);
end

dE_runtot = [0; cumsum(dE)];

% Hence find length of complete arc
dE_total = dE_runtot(end);

% Find N points equally spaced to span the supplied range
dE_spacing = dE_total./(n-1);

% Running total targets
dE_targets = dE_spacing * (0:(n-1));

% Linearly interpolate
lab = interp1(dE_runtot, lab, dE_targets, 'cubic');


switch lower(space)
    case 'lab'
        val = lab;
    case 'lch'
        val = lab2lch(lab);
    otherwise
        val = clab_colorspace([space '<-lab'],lab);
end

end