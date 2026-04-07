function Ktransp = ctranspose(K)
%
%    returns the conjugate transpose of the kronMatrix2 object K
%
%    K = A (x) B
%
%

A = K.a;
B = K.b;
Ktransp = kronMatrix2(A',B');