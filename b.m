% Evaluate the approximation of b at the point x given its truncated 
% decomposition on the Haar basis, Mu
function output = b(N,Mu,x)

H = zeros(N+2,2*N+1);
output = zeros(1,length(x));

for i = 1:length(x)
    for j = 1:N+2
        for m = 1:2*N+1
            H(j,m) = h(j-2,m-N-1,x(i));
        end
    end
    output(i) = sqrt(2) * Mu(1,:)*H(1,:)' + trace(Mu(2:N+2,:)*H(2:N+2,:)');
end

end