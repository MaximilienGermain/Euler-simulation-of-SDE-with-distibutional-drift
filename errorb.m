function err = errorb(Mu,N,Nmax)

diag(Mu(N+2:end,:)*Mu(N+2:end,:)')
err = sqrt(sum(diag(Mu(N+2:end,:)*Mu(N+2:end,:)')));

end