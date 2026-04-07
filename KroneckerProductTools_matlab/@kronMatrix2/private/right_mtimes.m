function y = right_mtimes(x, K)

%  y = right_mtimes(x, K)
%  Helper function for mtimes.m
%  K is a kronMatrix2 object, x is a numeric array (scalar,vector,or matrix)
%
%  Notice that y = (K.'*x.').'
%
y = (K.' * x.').';