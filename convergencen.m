function convergencen(min,max,X0,startNT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive)

frames = [];

for NT=min:max
    [~,~,frame,~,control,~] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive);
    frames = [frames frame];
end

video('Convergence time grid',frames)
close all

end