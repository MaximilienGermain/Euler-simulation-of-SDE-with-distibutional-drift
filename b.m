% Evaluate the approximation of b at the point x given its truncated 
% decomposition on the Haar basis, Mu
% Compute the approximation of b evaluated in a vector x,
% given its coefficients Mu on the N truncated Haar basis
function output = b(Mu,K,x)

sizeMu = size(Mu);
N = sizeMu(1) - 2;

% Construction of the evaluation of the haar functions at the vector x
hm = zeros(N+2,K*2^(N+1));
output = zeros(1,length(x));

for i = 1:length(x)
    for j = 1:N+2
        if (j==1)
            for m = 1:K*2
                hm(j,m) = h(j-2,m-K-1,x(i));
            end
        else
            for m = 1:K*2^(j-1)
                hm(j,m) = h(j-2,m-K*2^(j-2)-1,x(i));
            end
        end
    end
    hm = sparse(hm);
    % Calculus of the approximation of b(x)
    output(i) = sqrt(2) * Mu(1,:)*hm(1,:)' + trace(Mu(2:N+2,:)*hm(2:N+2,:)');
end

end