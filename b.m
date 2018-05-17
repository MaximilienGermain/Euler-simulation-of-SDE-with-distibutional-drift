% Evaluate the approximation of b at the point x given its truncated 
% decomposition on the Haar basis, Mu
% Compute the approximation of b evaluated in a vector x,
% given its coefficients Mu on the N truncated Haar basis
function output = b(Mu,K,x)

sizeMu = size(Mu);
N = sizeMu(1) - 2;

% Construction of the evaluation of the haar functions at the vector x
H = zeros(N+2,K*2^(N+1));
output = zeros(1,length(x));

for i = 1:length(x)
    for j = 1:N+2
        if (j==1)
            for m = 1:K*2
                H(j,m) = h(j-2,m-K-1,x(i));
            end
        else
            for m = 1:K*2^(j-1)
                H(j,m) = h(j-2,m-K*2^(j-2)-1,x(i));
            end
        end
    end
    H = sparse(H);
    % Calculus of the approximation of b(x)
    output(i) = sqrt(2) * Mu(1,:)*H(1,:)' + trace(Mu(2:N+2,:)*H(2:N+2,:)');
end

end