function Kout = uminus(Kin)
%
% Given a kronMatrix2 object, Kin, return Kout = -Kin
%
Kout = kronMatrix2(-Kin.a, Kin.b);