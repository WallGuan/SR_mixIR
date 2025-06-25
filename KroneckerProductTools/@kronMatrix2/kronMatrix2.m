function K = kronMatrix2( varargin )

%
% K = kronMatrix2( varargin )
%
% Constructor to create an Kronecker product matrix object. A more 
% extensive version of this, for sums of Kronecker products, can be
% found in IRtools.
%
% Object's fields are:
%    K.a = matrix A, where K = A (x) B
%    K.b = matrix B, where K = A (x) B
%
% Calling Syntax:
%       K = kronMatrix
%       K = kronMatrix(K);
%       K = kronMatrix(A, B);
%
%   where
%       K = a pre-existing kronMatrix object
%       A, B = both are matrices 
%

%  5/24/25 J. Nagy


%
% build the kronMatrix based on number and type of input arguments
%
switch nargin

  case 0
   K.a = [];
   K.b = [];
   K = class(K, 'kronMatrix2');
 
  case 1
   if isa(varargin{1}, 'kronMatrix2')
      K = varargin{1};
   else
      error('Wrong input argument')
   end

  case 2

   if isa(varargin{1}, 'double') & isa(varargin{2}, 'double')
     K.a = varargin{1};
     K.b = varargin{2};
     K = class(K, 'kronMatrix2');
   else
     error('Wrong input arguments')
   end % case 2

  otherwise
    error('Wrong input arguments')

end 