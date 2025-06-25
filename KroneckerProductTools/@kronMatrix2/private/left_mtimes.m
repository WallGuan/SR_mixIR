function y = left_mtimes(K, x)

%  y = left_mtimes(K, x)
%  Helper function for mtimes.m
%  K is a kronMatrix2 object, x is a numeric array (scalar,vector,or matrix)
%
%  if x is a scalar, use the property that
%  K * x = kron(A, B*x) 
%
%  if x is a vector such that x = vec(X),
%  use the property that 
%  kron(A, B) * x = vec(B * X * A')
%
%  if x has multiple columns, compute
%     y = [K*x(:,1) K*x(:,2) ... ]
%

if isscalar(x)
    y = kronMatrix2(K.a, x*K.b);
elseif size(x,1) == size(K.a,2)*size(K.b,2)
    % allocate space for y
    y = zeros(size(K.a,1)*size(K.b,1), size(x,2));
    % now sweep through columns of x 
    for j = 1:size(x,2)
        xj = reshape(x(:,j), size(K.b,2), size(K.a,2));
        y(:,j) = reshape(K.b * xj * K.a', size(K.a,1)*size(K.b,1), 1);
    end
else
    error('Dimension mismatch, incompatible array sizes')
end

