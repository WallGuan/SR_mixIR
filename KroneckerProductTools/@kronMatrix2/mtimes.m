function Kout = mtimes(K, M)
% N = mtimes(K, M)
%
%  kronMatrix2 multiplication;
%     multiply a kronMatrix2 by a matrix, a vector, or by
%     another kronMatrix2 (if possible),
%     

if (isa(K, 'kronMatrix2') && isnumeric(M))
    Kout = left_mtimes(K, M);
elseif (isnumeric(K) && isa(M, 'kronMatrix2'))
    Kout = right_mtimes(K, M);
elseif (isa(K, 'kronMatrix2') && isa(M, 'kronMatrix2'))
    if (size(K.a,2) == size(M.a,1) && size(K.b,2) == size(M.b,1))
        Anew = K.a * M.a;
        Bnew = K.b * M.b;
        Kout = kronMatrix2(Anew, Bnew);
    else
        error('Kron factors must be of compatible sizes for multiplication')
    end
else
    error('Wrong input arguments')
end

