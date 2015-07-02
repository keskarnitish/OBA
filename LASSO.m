classdef LASSO < handle
    properties
        XF
        X
        y
        msave
        num_features
        n
    end
    methods
        function obj = LASSO(X,y)
            obj.X = X;
            obj.y = y;
            obj.XF = X;
            obj.n = length(y);
            obj.num_features = size(X,2);
        end
        function f = func(obj)
            f = 0.5*norm(obj.msave-obj.y)^2;
        end
        
        function g = grad(obj)
            g = obj.X'*(obj.msave-obj.y);
        end
        
        function Hv = hess(obj,v)
            Hv = obj.XF'*(obj.XF*v);
        end
        
        function [] = setSave(obj,w)
            obj.msave = obj.X*w;
        end
        
        
        function [] = setMemoi(obj,F)
            obj.XF = obj.X(:,F);
        end
        
        function [] = removeMemoi(obj,to_delete_index)
            obj.XF(:,to_delete_index)=[];
        end
        
    end
end