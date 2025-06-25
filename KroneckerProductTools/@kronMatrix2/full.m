function A = full(K)
%
%                A = full(K);
%
%  Converts a kronMatrix2 matrix K to full storage organization.  
%
A = kron(K.a, K.b);
 