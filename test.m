X0 = 1;
NT = 1500;
N = 4;
K = 3;
T = 3;
%H = 0.95;
%H = 0.55;
H = 0.985;

%[random,xgrid,B,M] = fbm(H,K,N);
[xgrid,B,M] = createfBm(H,K,N);

X = eulerMethod(X0,NT,N,T,H,B,xgrid,0,K);