function [Hv] = reduced_HV_from_memory(V,PW,theta,invM,F,n)
    tV = V;
    V = zeros(n,1);
    V(F)= tV;
   
    %Hv = Hvfun(X,Y-X,save);
    if(~isempty(PW))
        PWV = PW'*V;
        Hv = theta*V - PW*(invM\PWV);
    else
        Hv = theta*V;
    end
    
    Hv = Hv(F);
end