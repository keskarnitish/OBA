classdef LogReg < handle
    properties
        XF
        b
        X
        y
        msave
        num_features
        n
    end
    methods
        function obj = LogReg(X,y)
            obj.X = X;
            obj.y = y;
            obj.XF = X;            
            obj.n = length(y);
            obj.num_features = size(X,2);
        end
        function f = func(obj)
            f = 1/obj.n*sum(log(1+obj.msave));
        end
        
        function g = grad(obj)
            g = -1/obj.n*((obj.b.*obj.y)'*obj.X)';
        end
        
        function Hv = hess(obj,v)
            c = obj.XF*v;
            C = c.*(obj.b-(obj.b).^2);
            Hv = 1/obj.n*(C'*obj.XF)'+1e-8*v;
        end
        
        function [] = setSave(obj,w)
            obj.msave = exp(-1*(obj.y).*(obj.X*w));
            obj.b = (obj.msave)./(1+obj.msave);
        end
        
        
        function [] = setMemoi(obj,F)
            obj.XF = obj.X(:,F);
        end
        
        function [] = removeMemoi(obj,to_delete_index)
            obj.XF(:,to_delete_index)=[];
        end
        
    end
end