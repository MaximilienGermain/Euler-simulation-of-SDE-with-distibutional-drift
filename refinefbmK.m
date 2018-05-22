function [newrandom,newgrid,newB,newM] = refinefbmK(random,xgrid,H,N,M,Ninit,previousB,p)

rng(ceil(100*N),'twister');
Nx = length(xgrid);
xgrid

%newgrid = linspace(min(xgrid)-1,max(xgrid)+1,Nx+2^(N+2))
newgrid = linspace(min(xgrid),max(xgrid)+1,Nx+2^(N+1))
min(xgrid)-1;
dx = 1/2^(N+1);

% Simulation of Nx - 1 independent gaussian variables
otherrandom1 = randn(1,2^(N+1))';
otherrandom2 = randn(1,2^(N+1))';

s = length(M)/2;
A = zeros(2*s,2^(N+1));
B = zeros(2^(N+1));

for i=1:2^(N+1)
    for j=1:2^(N+1)
        B(i,j) = dx^(2*H)*((i+2*s)^(2*H)+(j+2*s)^(2*H)-abs(i-j)^(2*H))/2;
    end
end
for i=1:(Ninit-1)
    for j=1:1:2^(N+1)
        A(i,j) = dx^(2*H)*((2^p*i)^(2*H)+(j+2*s)^(2*H)-abs(j+2*s-2^p*i)^(2*H))/2;
    end
    %dx*(2^p*i)
end
for k=1:p
    for i=2^(k-1)*(Ninit-1)+1:2^(k)*(Ninit-1)
        for j=1:2^(N+1)
            A(i,j) = dx^(2*H)*((2^(p-k)*(2*(i-2^(k-1)*(Ninit-1))-1))^(2*H)...
            +(j+2*s)^(2*H)-abs(j+2*s-2^(p-k)*(2*(i-2^(k-1)*(Ninit-1))-1))^(2*H))/2;            
        end
        k
        dx*(2^(p-k)*(2*(i-2^(k-1)*(Ninit-1))-1))
    end
end

C = A'/(M');
gamma = [M*M' A; A' B]
eig(gamma);
D = chol(B-C*C')';
M = [M zeros(2*s,2^(N+1)); C D];
newM = M;

newrandom = [random ; otherrandom1];
newvalues = M*newrandom;
%newvalues = M*other

newB=B;
newB=[previousB;newvalues(s+1:s+2^(N+1))];
newM=M;
newrandom=random;

end