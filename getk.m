% Index in the grid of size N + 1
function k = getk(j,m,N,K)

%k = 1 + m*2^(N+1-j) + N*2^(N+1);
k = 1 + 2^(N+1)*K + 2^(N+1-j)*m;

end