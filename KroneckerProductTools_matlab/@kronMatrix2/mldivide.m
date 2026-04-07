function x = mldivide(K,y)
%  kronMatrix2/mldivide
%
%     solves the linear system of equations
%     Kx = y 
%     when K = A (x) B
%

A = K.a;
B = K.b;
[ma, na] = size(A);
[mb, nb] = size(B);
Y = reshape(y, mb, ma);
Z = A \ Y';
X = B \ Z';
x = reshape(X, nb*na, 1);
