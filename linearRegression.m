function [beta0,beta1] = linearRegression(x,y)

n = length(x);
X = [ones(n,1) x];
beta = (X'*X)\(X'*y);
beta0 = beta(1);
beta1 = beta(2);

end