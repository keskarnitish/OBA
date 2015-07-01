function [Hv] = HV_from_memory(V,PW,theta,invM)

   
    %Hv = Hvfun(X,Y-X,save);
    if(~isempty(PW))
        PWV = PW'*V;
        Hv = theta*V - PW*(invM\PWV);
    else
        Hv = theta*V;
    end


end