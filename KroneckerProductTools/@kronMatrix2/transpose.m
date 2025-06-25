function Kout = transpose(Kin)

%  kronMatrix2 k.'
%    returns the transpose of the kronMatrix object k
%
%    Kin = A (x) B ==> Kin.' = A.' (x) B.'
%
%

Kout = kronMatrix2((Kin.a).',(Kin.b).');