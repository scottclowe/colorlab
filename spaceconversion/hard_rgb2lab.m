function varargout = hard_rgb2lab(varargin)
% Just needs to be a wrapper

[varargout{1:nargout}] = rgb2lab(varargin{:});

end