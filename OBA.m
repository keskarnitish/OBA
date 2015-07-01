function X = OBA(fun,lambda,options)

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


delta=1e3;


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
    
    fun.setXF(F);
    cnt=0;
    
    % Correction Strategy
    total_variables_moved=0;
    while(~isempty(Vk))
        cnt=cnt+1;
        if (quasi_newton == 1)
            if(mem.m>0)
                W = lbfgsFull(-PG,mem.DP(:,1:mem.m),mem.DG(:,1:mem.m),dot(mem.DG(:,mem.m),mem.DP(:,mem.m))/dot(mem.DG(:,mem.m),mem.DG(:,mem.m))   );
            else
                W = -PG;
            end
            W(setdiff((1:n)',F)) = 0;
            
            [PW,theta,invM] = lbfgs(mem.DP,mem.DG,mem.m,mem.updates,mem.size);
            %W = reduced_invHV_lbfgs(-PG,mem.DP,mem.DG,mem.m,mem.updates,mem.size,F);
            Hvfun = @(v) reduced_HV_from_memory(v,PW,theta,invM,F) ;
        else
            Hvfun = @(v) fun.hess(v);
            W = cg_steihaug ( Hvfun, PG, delta, A, CGrestol, maxCGiter);
        end
        
        
        XN=X+W;
        Vk = intersect(find(X==0 & XN.*Z<0),F);
        %HUGE_TABLE_GENERATOR
        A=union(A,Vk);
        
        to_delete_index = ismember(F,Vk);
        fun.removeXF(to_delete_index);
        
        %diag_of_H(to_delete_index)=[];
        F=setdiff(F,Vk);
        total_variables_moved=total_variables_moved+length(Vk);
    end
    
    % Projected Backtracking Line Search
    my_alpha = 1.0;
    q2 = @(XN,Hv) grad'*(XN-X)+(XN(F)-X(F))'*Hv*0.5+lambda*norm(XN,1);
    while(1)
        XN = X + my_alpha*W;
        %indices violating orthant bounds are projected to zero
        XN(F(Z(F).*XN(F)<0)) = 0;
        Hv_N = Hvfun(XN(F)-X(F));
        ls_lhs = -q2(XN,Hv_N) + lambda*norm(X,1);
        
        if(ls_lhs >= 0)
            X = XN;
            %f = q2(XN);
            %l1 = norm(XN,1);
            
            gradN = grad;
            gradN(F) = Hv_N + grad(F);
            grad = gradN;
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
        warning('OBA: ISTA Step Did Not Yield Descent. Increase Lipschitz Constant Estimate.');
    end
    
    X_trial = X_ista + alpha_dogleg*dogleg_d;
    fun.setSave(X_trial);
    
    f_trial = fun.func();
    
    while(1)
        if(f_trial+lambda*norm(X_trial(:),1)<=f_ista+lambda*l1_ista)
            grad_n = fun.grad();
            if(quasi_newton == 1)
                y = grad_n-grad;
                s = alpha_dogleg*dogleg_d;
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
                else
                    disp('Skipping');
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