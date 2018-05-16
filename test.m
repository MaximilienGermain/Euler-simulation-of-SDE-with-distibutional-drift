X0 = 1;
NT = 500;
N = 7;
T = 3;
H = 0.95;
[xgrid,B] = fbm(H,N);
Nx = 1 + (N+1)*2^(N+2);

X = eulerMethod(X0,NT,N,Nx,T,H,B,xgrid,0);