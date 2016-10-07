function x = mmean(x,dims,flag)

% -------------------------------------------------------------------------
% Input handling
if nargin<2
    % User specifies no dims and no flag -> max over everything
    dims = [];
    flag = 'x';
elseif nargin<3;
    % User specifies dims but no flag -> max over dims
    flag = '';
end

% -------------------------------------------------------------------------
if strcmpi(flag,'x');
    dims = setdiff(1:ndims(x),dims);
end

% -------------------------------------------------------------------------
for i=1:length(dims)
    x = mean(x,dims(i));
end

end