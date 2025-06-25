function varargout = svd(K)
%
% [U,S,V] = svd(K);
% s = svd(K);
%
% Computes an svd of the matrix represented of a kronMatrix2 object 
% %    K = A (x) B  
%
% On entry: K = a kronMatrix2 object
% On exit: U, V, S = kronMatrix2 objects, or 
%          s = vector of singular values 
%            ... such that U*S*V' is an approximation to G


[Ua, Sa, Va] = svd(K.a);
[Ub, Sb, Vb] = svd(K.b);

if (nargout == 3)
    varargout{1} = kronMatrix2(Ua, Ub);
    varargout{2} = kronMatrix2(Sa, Sb);
    varargout{3} = kronMatrix2(Va, Vb);
elseif (nargout == 1)
    s = kron(diag(Sa), diag(Sb));
    s = sort(s,'descend');
    varargout{1} = s;
else
    error('Incorrect number of output arguments.')
end

