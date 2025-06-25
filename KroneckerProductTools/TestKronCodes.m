%
% Some basic tests and illustrations of using the kronMatrix2 structure.
% 
% The kronMatrix2 structure provides a convenient way to store and do
% basic computations with a matrix K that is the Kronecker product of
% two matrices, K = A (x) B.
%
% The basic call K = kronMatrix2(A, B) builds the structure, with
% two fields: K.a = A, K.b = B. 
% Here's a quick example, where A and B are two random matrices:
A = rand(5,4);
B = rand(7,9);
K = kronMatrix2(A, B);
%
% There are a few basic operations coded for the kronMatrix2 object.
% One is full(K), which explicitly builds the Kronecker product into
% a full matrix.
%
% Note that MATLAB also has a function called kron that can explicitly
% build a Kronecker product. You don't want to use it if A and B are large,
% but it can be used to check things. For example:
%
K2 = kron(A, B);
norm(K2 - full((K)))
%
% Also notice that the kronMatrix2 object, K, uses a lot less storage than
% the full matrix K2:
%
whos K K2
%
% Another important operation coded for the kronMatrix2 object is 
% matrix multiplciation (mtimes). Here's a simple example to use it,
% and also to test correctness:
%
x = rand(size(K,2),1);
b = K*x;
norm(b - K2*x)
%
% It should also be possible to multiply K to a matrix x with multiple
% columns:
%
x = rand(size(K,2),10);
b = K*x;
norm(b - K2*x)
%
% It should also be possible to multiple on the left with a row vector,
% or with a matrix of compatible dimensions:
%
y = rand(10,size(K,1));
z = y*K;
norm(z - y*K2)
%
% Solving linear systems with backslash can also be done pretty easily.
% But for this test, we shoudl build a square, nonsingular matrix.
%
A = rand(4,4); 
B = rand(7,7);
K = kronMatrix2(A, B);
K2 = kron(A, B);
b = rand(size(K,1),1);
x = K\b;
norm(x - K2\b)/norm(x)
%
% Transpose should also work:
%
norm(full(K') - K2')
%
% Finally, it's nice to have the SVD working for a kronMatrix2 object.
% You can compute the SVD using:
%
[U, S, V] = svd(K);
%
% But notice that U, S, V are kronMatrix2 objects. You don't want to 
% explicitly build these for large problems, but we can build them
% for smaller problems. 
%
whos U S V
% 
% We can also just compute a vector with the singular values:
%
s = svd(K);
figure(1), clf
axes('FontSize', 18), hold on
plot(s, 'b-o', 'LineWidth', 2)
title('Singular Values of K', 'FontSize', 18)