[K, g, f] = heat(256);
[U,S,V] = svd(K);
s = diag(S);
subplot(2,2,1), plot(s, 'o')
subplot(2,2,2), plot(V(:,1))
subplot(2,2,3), plot(V(:,20))
subplot(2,2,4), plot(V(:,100))