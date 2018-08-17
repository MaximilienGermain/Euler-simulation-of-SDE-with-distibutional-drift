function err = errorb(Muerr,N,Nmax)

diag(Muerr(N+2:end,:)*Muerr(N+2:end,:)')
err = sqrt(sum(diag(Muerr(N+2:end,:)*Muerr(N+2:end,:)')));

end