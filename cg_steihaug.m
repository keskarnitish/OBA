function [ p, num_cg, iflag ] = cg_steihaug ( Hv, b, delta, A, errtol, maxit)
% CG-Steihaug based on 
% http://users.eecs.northwestern.edu/~morales/PSfiles/cg_steihaug.m
n      = length(b);  
%errtol = params(1); maxit  = params(2);
iprnt  = 0;
iflag  = ' ';
%
g  = b; 
g(A) = 0;
x  = zeros(n,1); r = -g;
%
z   = r;
rho = z'*r;   
tst = norm(r,'inf'); 
flag  = ''; 
terminate = errtol*norm(g,'inf');   it = 1;    hatdel = delta*(1-1.d-6);
rhoold = 1.0d0;
if iprnt > 0 
    fprintf(1,'\n\tThis is an output of the CG-Steihaug method. \n\tDelta = %7.1e \n', delta);
    fprintf(1,'   ---------------------------------------------\n');
end
flag = 'We do not know ';
if tst <= terminate; flag  = 'Small ||g||    '; end 
    F = (1:n)';
    F(A) = [];
    w = x;
while((tst > terminate) & (it <=  maxit) & norm(x) <=  hatdel)
    %
    if(it == 1) 
        p = z;
    else
        beta = rho/rhoold;
        p = z + beta*p;
    end
    %
    %
    w = zeros(size(g));
    w(F)  = Hv(p(F));
    
    alpha = w'*p;
    %
    % If alpha < = 0 head to the TR boundary and return
    %
    ineg = 0;
    if(alpha <=  0)
        ac = p'*p; bc = 2*(x'*p); cc = x'*x - delta*delta;
        alpha = (-bc + sqrt(bc*bc - 4*ac*cc))/(2*ac);
        flag  = 'negative curvature';
        iflag = 'NC';
    else
        alpha = rho/alpha;
        if norm(x+alpha*p) > delta
            ac = p'*p; bc = 2*(x'*p); cc = x'*x - delta*delta;
            alpha = (-bc + sqrt(bc*bc - 4*ac*cc))/(2*ac);
            flag  = 'boundary was hit';
            iflag = 'TR';
        end
    end
    x   =  x + alpha*p;
    r   =  r - alpha*w;
    tst = norm(r,'inf');
    if tst <= terminate; flag = '||r|| < test   '; iflag = 'RS'; end;
    if norm(x) >=  hatdel; flag = 'close to the boundary'; iflag = 'TR'; end
    
    if iprnt > 0 
        fprintf(1,' %3i    %14.8e   %s  \n', it, tst, flag);
    end
    rhoold = rho;
    z   = r; 
    rho = z'*r;
    it  = it + 1;
end %           
if it > maxit; iflag = 'MX'; end;

num_cg = it;
p = x;
%