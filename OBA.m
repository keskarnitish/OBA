function X = OBA(fun,lambda,options)
% OBA : A Second-Order Method for Convex L1-Regularized Optimization with Active Set Prediction
%  Minimization algorithm intended for composite convex functions and
%   specifically for high-dimensional problems.
%  Usage: X = OBA(fun,lambda,[options])
%  returns the approximate minimizer of the function fun + lambda*l1
%  where fun is object with member functions for function, gradient and
%  Hessian (through Hessian-vector products) evaluations.
%
%   Input parameters
%    fun is a required object, with 3 required member functions (for L-BFGS) and 5
%    (for Newton's method). For ease of design, please use the funTemplate
%    class to design your own class.
%    Logistic Regression and LASSO are included in this package.
%
%    If using quasi-Newton, the following three functions are required.
%       fun.setSave(X): Given an iterate X, this function sets X as the
%       iterate for the object. Any check-pointing computations can be stored
%       in the object.
%
%       fun.func(): Having saved an iterate, this function computes the
%       function value at the saved point.
%
%       fun.grad(): Having saved an iterate, this function computes the
%         gradient value at the saved point.
%
%    For Newton, in addition to the above functions, the following two
%    member functions are also required.
%
%       fun.setMemoi(F) and removeMemoi(del): Since OBA solves QPs in a subspace,
%       it can benefit from memoization with enables subspace-Hessian
%       computations. setMemoi(F) uses the free-set F to create check-pointed
%       computations which enable easier subspace-Hessian computations and
%       removeMemoi(del) removes the effect of the elements of del from the
%       check-pointed computations. Please refer to LASSO.m for example.
%
%       fun.hess(v): Given a (subspace) vector v and having saved an iterate, this
%       function computes the subspace Hessian-vector product at the saved point.
%
%    options is an optional struct. In order to specify options, the user
%    should first run options = GenOptions() and modify the structure as
%    necessary. The details about the structure and the default values are
%    as follows.
%
%      options.optol: termination tolerance
%          (default: 1e-6)
%      options.qn: Quasi-Newton, 0 (Newton's Method), or 1 (quasi-Newton)
%          (default: 0)
%      options.mem_size: quasi-Newton memory size
%          (default: 20)
%      options.maxiter: max number of iterations
%          (default: 1000)
%      options.printlev: print level, 0 (no printing) or 1
%          (default: 1)
%      options.CGtol: CG termination tolerance (for Newton's Method)
%          (default: 1e-1)
%      options.maxCGiter: max number of CG iterations (Newton's Method)
%          (default: 1000).
%   References:
%   A Second-Order Method for Convex L1-Regularized Optimization with
%   Active Set Prediction. Link: http://arxiv.org/abs/1505.04315
%
%   Send comments/bug reports to Nitish Shirish Keskar,
%   keskar.nitish@u.northwestern.edu
%   Version 1.1, 2015, see GPL license info below.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  OBA 1.1 Copyright (C) 2015 Nitish Shirish Keskar
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization

if nargin<3
    options = GenOptions();
end

maxit        = options.maxiter;
termtol      = options.optol;
printlevel   = options.printlev;
CGrestol     = options.CGtol;
maxCGiter    = options.maxCGiter;
quasi_newton = options.qn;
n = fun.num_features;
X = zeros(n,1);
perc=min(99,100-1/length(X)*100);



% Quasi Newton

if(quasi_newton==1) %initialize memory for LBFGS
    %H = eye(n);
    mem_size = options.mem_size;
    mem = []; mem.updates=0;  mem.size=0;
    mem.DP = zeros(n,mem.size);
    mem.DG = zeros(n,mem.size);
    mem.m = 0;
    mem.updates=0;
    mem.size = mem_size;
end


%% FUCTION, L1 AND GRADIENT EVALUATION AT INITIAL POINT
fun.setSave(X);
l1 = norm(X(:),1);
f=fun.func();
grad=fun.grad();


%% Outer (Optimization) Loop
iteration = 0;
while(1)
    
    % Check Optimality Conditions
    tt = 1e-8;
    X_ista = soft_threshold_l1(X-tt*grad,lambda*tt);
    opterr = norm(max(min(grad+lambda*ones(size(X)),X),grad-lambda*ones(size(X))),'inf');
    
    if(printlevel == 1)
        fprintf('OBA: iter: %3d\t opt_err: %e f: %3.4e |X|1: %7.4f obj: %3.4e |X|0: %d \n',...
            iteration , opterr , full(f) , full(l1), ...
            full(f+lambda*l1), full(sum(sum(X~=0))));
    end
    
    
    if(opterr<termtol)
        disp('OBA: Optimal Solution Found.');
        break;
    end
    
    iteration = iteration + 1;
    if(iteration>maxit)
        disp('OBA: Maximum Number of Iterations Exceeded.');
        break;
    end
    
    % Preliminary Computations
    
    %orthant indicator
    Z=SteepestDescent_ZeroActive(X,grad,lambda);
    
    %pseudo-gradient (Minimum Norm Subgradient)
    PG = grad + lambda*Z;
    Vk=0;
    A = find( X==0 & abs(grad)<lambda ); % A^k
    F_zeros = find(X==0 & abs(grad)>=lambda); % U^k
    F_nonzeros = (1:n)';
    F_nonzeros(union(A,F_zeros))=[]; % F^k
    [~, IND_OF_F] = sort(abs(PG(F_zeros))); % Sort U^k by PG
    F_zeros=F_zeros(IND_OF_F);
    
    % Selection Strategy
    tau_hat = floor((1-perc/100)* (length(X)));
    tau = min(tau_hat,length(F_zeros));
    F = union (F_nonzeros, F_zeros(end-tau+1:end));
    F_zeros(end-tau+1:end)=[];
    A=union(A,F_zeros);
    
    fun.setMemoi(F);
    cnt=0;
    
    % Correction Strategy
    total_variables_moved=0;
    while(~isempty(Vk))
        if (quasi_newton == 1)
            [PW,theta,invM] = subspaceSolvers.lbfgs(mem.DP,mem.DG,mem.m,mem.updates,mem.size);
            Hvfun = @(v) subspaceSolvers.reduced_HV_from_memory(v,PW,theta,invM,F,n) ;
            W = subspaceSolvers.reduced_invHV_from_memory(-PG,PW,theta,invM,F);
        else
            Hvfun = @(v) fun.hess(v);
            W = subspaceSolvers.CG ( Hvfun, PG, A, CGrestol, maxCGiter);
        end
        
        
        XN=X+W;
        Vk = intersect(find(X==0 & XN.*Z<0),F);
        %HUGE_TABLE_GENERATOR
        A=union(A,Vk);
        
        to_delete_index = ismember(F,Vk);
        fun.removeMemoi(to_delete_index);
        
        %diag_of_H(to_delete_index)=[];
        F=setdiff(F,Vk);
        total_variables_moved=total_variables_moved+length(Vk);
    end
    
    % Projected Backtracking Line Search
    X_old = X;
    my_alpha = 1.0;
    q2 = @(XN,Hv) grad'*(XN-X)+(XN(F)-X(F))'*Hv*0.5+lambda*norm(XN,1);
    while(1)
        XN = X + my_alpha*W;
        %indices violating orthant bounds are projected to zero
        XN(F(Z(F).*XN(F)<0)) = 0;
        %fun.setSave(XN);
        Hv_N = Hvfun(XN(F)-X(F));
        ls_lhs = -q2(XN,Hv_N) + lambda*norm(X,1);
        
        if(ls_lhs >= 0)
        %if( (f+lambda*l1) > (fun.func() + lambda*norm(XN,1)))
            X = XN;
            %f = q2(XN);
            %l1 = norm(XN,1);
            
            %             gradN = grad;
            %             gradN(F) = Hv_N + grad(F);
            %             grad = gradN;
            break;
            
        else
            my_alpha = my_alpha/2;
        end
    end
    
    % Adapt Selection Parameter
    if(total_variables_moved==0)
        tmp_perc = 100-perc;
        tmp_perc = tmp_perc * 2;
        perc = 100-tmp_perc;
        perc = max(0,perc);
        %fprintf('Decreasing Perc to %f \n',perc);
    end
    
    
    % ISTA Globalization Strategy
    dogleg_d = X-X_ista;
    alpha_dogleg = 1;
    l1_ista = norm(X_ista(:),1);
    f_ista = f+grad'*(X_ista-X)+1/(2*tt)*norm(X_ista-X)^2;
    
    if(f+lambda*l1-f_ista-lambda*l1_ista>=-1e-8)
        %warning('OBA: ISTA Step Did Not Yield Descent. Increase Lipschitz Constant Estimate.');
    end
    
    X_trial = X_ista + alpha_dogleg*dogleg_d;
    fun.setSave(X_trial);
    
    f_trial = fun.func();
    
    while(1)
        if(f_trial+lambda*norm(X_trial(:),1)<=f_ista+lambda*l1_ista)
            grad_n = fun.grad();
            
            if(quasi_newton == 1)
                y = grad_n+lambda*SteepestDescent_ZeroActive(X_trial,grad_n,lambda)-PG;
                s = X_trial - X_old;
                if (dot(y,s)>1e-8)
                    %rho = 1/dot(s,y);
                    %H = (eye(n) - rho*s*y')*H*(eye(n) - rho*y*s') + rho*(s*s');
                    if mem.updates < mem.size
                        mem.DG(:,mem.updates+1) = y;
                        mem.DP(:,mem.updates+1) = s;
                    else
                        mem.DG = [mem.DG(:,2:mem.size) y];
                        mem.DP = [mem.DP(:,2:mem.size) s];
                    end
                    mem.updates=mem.updates+1;
                    mem.m = min(mem.size,mem.updates);
                    %else
                    %disp('Skipping');
                end
            end
            X = X_trial;
            f = f_trial;
            grad = grad_n;
            l1 = norm(X(:),1);
            break;
        else
            alpha_dogleg = alpha_dogleg*0.5;
            if(alpha_dogleg<1e-4)
                X = X_ista ;
                l1 = norm(X_ista,1);
                fun.setSave(X_ista);
                f_ista=fun.func();
                f = f_ista;
                grad = fun.grad();
                break
            end
        end
        X_trial = X_ista + alpha_dogleg*dogleg_d;
        fun.setSave(X_trial);
        f_trial = fun.func();
        
    end
    
    
    
    
    
end

end

function output=soft_threshold_l1(y,lambda)
output=sign(y).*max(abs(y)-lambda,0);
end

function Z=SteepestDescent_ZeroActive(X,grad,lambda)
Z = zeros(size(grad));
Z(X>0 | (X==0 & grad<-lambda)) = 1; 
Z(X<0 | (X==0 & grad>lambda)) = -1;
end