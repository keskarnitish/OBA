function Z=SteepestDescent_ZeroActive(X,grad,lambda)
% OBM_COR
Z = zeros(size(grad));
Z(X>0 | (X==0 & grad<-lambda)) = 1;  %FGN:set conditions with a tolerance?
Z(X<0 | (X==0 & grad>lambda)) = -1;
end