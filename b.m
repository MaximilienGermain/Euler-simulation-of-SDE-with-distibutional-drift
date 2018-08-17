% Compute the approximation of b evaluated in a vector x,
% given its coefficients Mu on the N truncated Haar basis
function output = b(Mu,K,x)

sizeMu = size(Mu);
N = sizeMu(1) - 2;

% Construction of the evaluation of the haar functions at the vector x
output = zeros(1,length(x));

for i = 1:length(x)
    hm = zeros(N+2,K*2^(N+1));
    for j = 1:N+2
        if (j==1)
            if (floor(x(i)) < K) && (floor(x(i)) >= -K)
                hm(j,floor(x(i))+K+1) = h(j-2,floor(x(i)),x(i));
               % l(j,floor(x(i))+K+1) = h(j-2,floor(x(i)),x(i));
            end
        else
            p = floor(2^(j-2)*x(i));
            if (p < K*2^(j-2)) && (p >= -K*2^(j-2))
                 hm(j,p+K*2^(j-2)+1) = 2^(j-2)* h(j-2,p,x(i));
                 % hm(j,p+K*2^(j-2)+1) = h(j-2,p,x(i)); ok
            end
        end
    end
    hm = sparse(hm);
    % Calculus of the approximation of b(x)
    output(i) = Mu(1,:)*hm(1,:)' + trace(Mu(2:N+2,:)*hm(2:N+2,:)');
end

end