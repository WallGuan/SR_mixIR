function varargout = size(K, dim)
%  kronMatrix2/size
%
%     returns the size of a kronMatrix2 object.
%

s = size(K.a) .* size(K.b);

if nargin == 2
    s = s(dim);
end
switch (nargout)
  case 0
    disp(sprintf('\nans =\n'));
    disp(s);
  case 1
    varargout{1} = s;
  case 2
    if nargin == 2
        error('Too many output arguments')
    else
        varargout{1} = s(1);
        varargout{2} = s(2);
    end
end
