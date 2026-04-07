function [A,b,x] = heat(n,kappa)
% Set default kappa.
if (nargin==1), kappa = 1; end
% Initialization.
h = 1/n; t = h/2:h:1;
c = h/(2*kappa*sqrt(pi));
d = 1/(4*kappa^2);
% Compute the matrix A.
k = c*t.^(-1.5).*exp(-d./t);
r = zeros(1,length(t)); r(1) = k(1); A = toeplitz(k,r);
% Compute the vectors x and b.
if (nargout>1)
  x = zeros(n,1);
  for i=1:n/2
    ti = i*20/n;
    if (ti < 2)
      x(i) = 0.75*ti^2/4;
    elseif (ti < 3)
      x(i) = 0.75 + (ti-2)*(3-ti);
    else
      x(i) = 0.75*exp(-(ti-3)*2);
    end
  end
  x(n/2+1:n) = zeros(1,n/2);
  b = A*x;
end