X0 = 1;
NT = 500;
N = 15;
T = 3;

Mu = zeros(N+2,2*N+1);
for m=1:2*N+1
    Mu(1,m) = (2*(m-N-1)+1)/(2*sqrt(2));
end    
    
for j=2:N+2
    Mu(j,:) = - 1/(2^(j));
end

X = eulerMethod(X0,NT,N,T,Mu);