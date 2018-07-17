function err = errorb(Mu,N)

diag(Mu(N+1:end,:)*Mu(N+1:end,:)')
err = sqrt(sum(diag(Mu(N+1:end,:)*Mu(N+1:end,:)')));

end