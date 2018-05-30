% Here newM is the the transpose of the cholesky root
function [newrandom,newgrid,newB,newM] = refinefbm(random,xgrid,H,N,M,Ninit,p,seed)

rng(seed*(p+1),'twister');
Nx = length(xgrid);
newgrid = linspace(min(xgrid),max(xgrid),2*Nx-1);
dx = 1/2^(N+1);

% Simulation of Nx - 1 independent gaussian variables
otherrandom = randn(1,Nx-1)';

s = length(M);
A = zeros(s);
B = zeros(s);

for i=1:s
    for j=1:s
        B(i,j) = dx^(2*H)*((2*i-1)^(2*H)+(2*j-1)^(2*H)-4^(H)*abs(i-j)^(2*H))/2;
    end
end
for i=1:(Ninit-1)
    for j=1:s
        A(i,j) = dx^(2*H)*((2^p*i)^(2*H)+(2*j-1)^(2*H)-abs(2*j-1-2^p*i)^(2*H))/2;
    end  
end
for k=1:p-1
    for i=2^(k-1)*(Ninit-1)+1:2^(k)*(Ninit-1)
        for j=1:s
            A(i,j) = dx^(2*H)*((2^(p-k)*(2*(i-2^(k-1)*(Ninit-1))-1))^(2*H)...
            +(2*j-1)^(2*H)-abs(2*j-1-2^(p-k)*(2*(i-2^(k-1)*(Ninit-1))-1))^(2*H))/2;
            %dx*(2^(p-k)*(2*(i-2^(k-1)*(Ninit-1))-1))
        end
    end
end
%M*M';
%A;
%B;
%C = A'/(M');
gamma = [M*M' A; A' B];
%eig(gamma);
M = chol(gamma)';
%D = chol(B-C*C')';
%M = [M zeros(s); C D];
newM = M;

newrandom = [random ; otherrandom];
newvalues = M*newrandom;

index = 2*(Ninit-1);
temp = fusion(newvalues(Ninit:index),newvalues(1:Ninit-1));
while index < length(newvalues)
    temp = fusion(newvalues(index+1:2*(index)),temp);
    index=2*index;
end
newB=[0 ;temp];