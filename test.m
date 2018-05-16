X0 = 0.2;
NT = 50;
N = 2;
T = 1;

Mu = zeros(N+2,2*N+1);
for m=1:2*N+1
    Mu(1,m) = (2*(m-N-1)+1)/(2*sqrt(2));
end    
    
for j=2:N+2
    Mu(j,:) = - 1/(2^(j));
end

X = eulerMethod(X0,NT,N,T,Mu);