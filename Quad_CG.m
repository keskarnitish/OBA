function [X,evals,err,res] = Quad_CG(P,P_0,Hv,save,PG,A,tol,maxiter,precond)
%CG algorithm for computing an approximate Newton step
%Let M=C'*C be the preconditioner, where C is a vector standing for the
%diagonal entries of a diagonal matrix
evals = 0;
if(precond>0)
    C = diag(P); %M(A)=[];
end
R = -PG;
if(precond>0)
    R = repmat(C,1,length(C)).*R;
end
R(A)=0;   Q = R;
if(precond>0)
    Q = repmat(C,1,length(C)).*Q;
end


X = zeros(size(P));


cg_iteration = 0;
err=3;
%fprintf('\n IN CG MAXITER=%d\n',maxiter);
while(cg_iteration < maxiter)
    res = norm(R,'inf')/norm(PG,'inf');
    if(res<tol)
        err=0;
        break;
    end
    cg_iteration = cg_iteration + 1;

    Y = Hv(P,Q,save);
    evals = evals +1;
       
    %end
    Y(A) = 0;
    RR = R(:)'*R(:);
    alpha = RR / (Q(:)'*Y(:));
    X = X + alpha*Q;
    %norm(reshape(X,100,100))
    if(precond>0)
        R = R - alpha*repmat(C,1,length(C)).*Y;
    else
        R = R - alpha*Y;
    end
    beta =  (R(:)'*R(:)) / RR;
    if(precond>0)
        Q = repmat(C,1,length(C)).*R + beta*Q;
    else
        Q = R + beta*Q;
    end
end
%fprintf('CG:num_cg_iter:%d\t res:%e \n',cg_iteration,res);
end