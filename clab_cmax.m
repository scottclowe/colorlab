% CLAB_CMAX Maximises chroma for a set of colours
%    CLAB_CMAX(L,H) where L and H are vectors of the same length finds the 
%    maximum chroma which is in gamut for all the L(i),H(i) pairs with L as
%    Lightness and H as Hue in the CIELCH colorspace.
%     
%    CLAB_CMAX(L,H) where one of L and H is singleton finds the maximum
%    chroma for all L,H(i) or L(i),H pairs.
%    
%    CLAB_CMAX([],H) finds the maximum chroma for each hue in H under the
%    assumption that the colours can be specified at any lightness so long
%    as the lightness is the same for all colours.
%    
%    CLAB_CMAX(LCH) where LCH is a matrix of colours in CIE L*C*h* finds
%
%    CLAB_CMAX(SPACENAME,X) where X is a matrix of colours in the
%    colorspace SPACENAME. SPACENAME can be 'rgb', 'lab' or 'lch'.
%    Alternatively, SPACENAME can be 'lh', with X a matrix with Lightness
%    in the first column and Hue in the second (Chroma is omitted).
%    
%    CLAB_CMAX(...,GAMUT) maximises chroma within the provieded GAMUT
%    structure. If this is omitted, the default gamut is used as found by
%    FETCH_CIELCHAB_GAMUT.
%    
%    [X,C] = CLAB_CMAX(...) returns the common chroma, C, as well as the
%    chroma-maximised colours, X, in the same colorspace as the colorspace
%    of the input.

% Scott Lowe, 2014-03-05

% Notes:
%   UPLab is currently unsupported

function [X,C] = clab_cmax(varargin)

% Input handling ==========================================================
error(nargchk(1, 3, nargin, 'struct'));

nargin2 = nargin; % Number of inputs excluding gamut structure, if present
if isstruct(varargin{nargin})
    if nargin==1; error('Incorrect input type: struct'); end;
    gamut = varargin{nargin};
    nargin2 = nargin2-1;
else
    gamut = fetch_cielchab_gamut;
end

inspace  = 'lch';
dynamicL = false;

if nargin2==1
    lch = varargin{1};
    if size(lch,2)~=3; error('LCH input must be n-by-3'); end;
    
elseif ischar(varargin{1})
    if ~isnumeric(varargin{2})
        error('Second input must be numeric');
    end
    switch lower(varargin{1})
        case 'lh' % omitting chroma
            if size(varargin{2},2)~=2; error('Input must be n-by-2'); end
            lch = ones(size(varargin{2},1),3);
            lch(:,[1 3]) = varargin{2};
            inspace  = 'lch';
        case 'lch'
            if size(varargin{2},2)~=3; error('Input must be n-by-3'); end
            lch = varargin{2};
            inspace  = 'lch';
        case 'lab'
            if size(varargin{2},2)~=3; error('Input must be n-by-3'); end
            lch = lab2lch(varargin{2});
            inspace = 'lab';
        case 'rgb'
            if size(varargin{2},2)~=3; error('Input must be n-by-3'); end
            lch = lab2lch(rgb2lab(varargin{2}));
            inspace = 'rgb';
        otherwise
            error('Unfamiliar colorspace %s',varargin{1});
    end
    
elseif ~isnumeric(varargin{1})
    error('Incorrect input #1 type: %s',class(varargin{1}));
    
elseif ~isnumeric(varargin{2})
    error('Incorrect input #2 type: %s',class(varargin{2}));
    
elseif isempty(varargin{1})
    hvec = varargin{2};
    dynamicL = true;
    
elseif isequal(size(varargin{1}),size(varargin{2})) || ...
        isequal(size(varargin{1}),size(varargin{2}'))
    lch = zeros(numel(varargin{1}),3);
    lch(:,1) = varargin{1}(:); % Vector L
    lch(:,2) = 1;              % Ensure C>0. Will override with maxC.
    lch(:,3) = varargin{2}(:); % Vector H
    
elseif all(size(varargin{1})==1)
    lch = zeros(numel(varargin{2}),3);
    lch(:,1) = varargin{1};    % Scalar L
    lch(:,2) = 1;              % Ensure C>0. Will override with maxC.
    lch(:,3) = varargin{2}(:); % Vector H
    
elseif all(size(varargin{2})==1)
    lch = zeros(numel(varargin{1}),3);
    lch(:,1) = varargin{1}(:); % Vector L
    lch(:,2) = 1;              % Ensure C>0. Will override with maxC.
    lch(:,3) = varargin{2};    % Scalar H
    
else
    error('Unsure how to handle to inputs with sizes %s and %s',...
        mat2str(size(varargin{1})),mat2str(size(varargin{2})));
    
end

% Main script =============================================================
if dynamicL
    % Round to nearest h-interval in loaded gamut
    hvec = hvec(:);
    hvec = round(hvec/gamut.hintv)*gamut.hintv;
    C = gamut.lchmesh.cgrid(ismember(gamut.lchmesh.hvec,hvec),:);
    C = min(C,[],1);
    [C,I] = max(C,[],2);
    L = gamut.lchmesh.Lvec(I);
    
    lch = zeros(numel(hvec),3);
    lch(:,1) = L;
    lch(:,2) = C;
    lch(:,3) = hvec(:);
    
else
    % We get isingamut to do all the donkey work for us
    % It interpolates on the L-h grid to find the maximum C in gamut at each
    % point.
    % We do this for all the input colours
    [TF,lch] = isingamut(lch,gamut,'lch');
    
    % Then we find their minimum to get the maximum chroma which is in gamut at
    % all L-h locations
    C = min(lch(:,2));
    
    % Now set the chroma of each of our colors to be the maximum useable chroma
    lch(:,2) = C;
end

% Convert to output colorspace --------------------------------------------
switch inspace
    case 'lch'
        X = lch;
    case 'lab'
        X = lch2lab(lch);
    case 'rgb'
        X = hard_lab2rgb(lch2lab(lch));
    otherwise
        error('Unfamiliar output space');
end

end