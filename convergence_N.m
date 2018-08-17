function convergence_N(min,max,X0,startNT,startN,NT,T,H,testId,Kmax,seed,Nx,graphHaar,control,PlotActive)

frames = [];

for N=min:max
    [xgrid,B,~] = createfBm(H,Kmax,N,startN,Nx,-Kmax,1000);
    [Mu,~] = computeMu(B,N,testId,Kmax,H);
    [~,~,frame,~,~,~] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive);
    frames = [frames frame];
end
video('Convergence x grid',frames)
close all

end