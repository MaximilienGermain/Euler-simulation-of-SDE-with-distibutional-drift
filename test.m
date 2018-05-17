X0 = 1;
NT = 1500;
N = 2;
K = 5;
T = 3;
H = 0.95;
Nx = 1+K*2^(N+2);
[xgrid,B] = fbm(H,K,Nx,N);

X = eulerMethod(X0,NT,N,Nx,T,H,B,xgrid,0,K);