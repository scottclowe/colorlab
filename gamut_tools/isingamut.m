% Checks whether points are within the specified gamut
function [TF,P2] = isingamut(P,g,space)

% Parse inputs
switch lower(space)
    case 'lch'
        % Convert h so it is between [0,360)
        P(:,3) = mod(P(:,3),360);
        
    case 'lab'
        % Need to convert to Lch
        L = P(:,1);
        a = P(:,2);
        b = P(:,3);
        
        c = sqrt(a.^2+b.^2);
        h = atan2(b,a)/pi*180;
        h = mod(h,360);
        
        P = [L c h];
        
    case 'rgb'
        % Trivial case
        if size(P,2)==3
            dim = 2;
        elseif size(P,1)==3
            dim = 2;
        else
            error('bad input');
        end
        TF = any(P<0,dim) | any(P>1,dim);
        return;
        
    case {'xyy','xyl'}
        if size(P,2)==2 || size(P,2)==3
            dim = 2;
        else
            error('bad input');
        end
        Y = find_gamut_xy_maxY(P(:,1), P(:,2), g.space);
        
        TF = ~isnan(Y);
        if size(P,dim)==3;
            TF = TF | P(:,3)>Y;
        end
        
        P2 = P;
        P2(:,3) = Y;
        return;
    
    otherwise
        error('Unknown color space');
        
end

% Get a mesh version of the gamut
if ~isfield(g,'lchmesh')
    g.lchmesh = make_gamut_mesh(g);
end

% Add a 360 hue to this
if ~ismember(360,g.lchmesh.hvec) && ismember(0,g.lchmesh.hvec)
    g.lchmesh.Lgrid(end+1,:) = g.lchmesh.Lgrid(1,:);
    g.lchmesh.hgrid(end+1,:) = 360;
    g.lchmesh.cgrid(end+1,:) = g.lchmesh.cgrid(1,:);
end

% We can then linearly interpolate on this fine mesh
Cedge = interp2(g.lchmesh.Lgrid, g.lchmesh.hgrid, g.lchmesh.cgrid, P(:,1), P(:,3));

% Compare the input with the interpolated gamut chroma to see if it is in
TF = P(:,2) <= Cedge;

if nargout<2; return; end;

% We can also return the interpolated gamut point
P2 = [P(:,1) Cedge P(:,3)];

% This may need to be converted from Lch to a different space
switch lower(space)
    case 'lch'
        % Do nothing
        
    case 'lab'
        % Convert back to L*a*b*
        L2 = P2(:,1);
        c2 = P2(:,2);
        h2 = P2(:,3);
        
        a2 = c2 .* cos(h2/180*pi);
        b2 = c2 .* sin(h2/180*pi);
        
        P2 = [L2 a2 b2];
        
    otherwise
        error('Unknown color space');
end

end