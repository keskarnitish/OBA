classdef funTemplate < handle
    properties
        %Fill this with the datasets necessary to make function and
        %gradient evaluations. One possible "gotcha" is that the function
        %and gradient evaluations do NOT use the argument 
    end
    methods
        function obj = funTemplate(args)
            %Class constructor. Use this to define the properties listed
            %above. 
        end
        function f = func(obj)
            %Compute the function using the saved information (could be
            %iterate information or some check-pointed computation).
        end
        
        function g = grad(obj)
            %This function is used to define the gradient computation which
            %uses a syntax similar to the function evaluation.
        end
        
        function Hv = hess(obj,v)
            %This function computes a Hessian-vector product. This
            %function is optional and can be left empty if the user is only
            %interested in the L-BFGS algorithm. 
        end
        
        function [] = setSave(obj,w)
            %This function is required. Given the iterate w, this function
            %can be used to store the iterate information in the object or
            %to make some check-pointed computation which is used for the
            %various (function, gradient and Hessian) computations. 
        end
        
        
        function [] = setMemoi(obj,F)
            %Since OBA reduces the search space and then solves a (smooth)
            %QP, it can benefit from having check-pointed information which
            %allows easier Hessian-vector computations in the subspace.
            %This function, with the free set F provided, calculates an intermediate
            %(i.e. check-pointed) value which can be used for subspace
            %Hessian-vector products. 
        end
        
        function [] = removeMemoi(obj,to_delete_index)
            %This function is closely tied to the setMemoi function. The
            %setMemoi function calculates an intermediate value which can
            %be used for subspace Hessian-vector computation. As a part of
            %the correction mechanism, elements are removed from the free
            %set. This function uses the computations made by setMemoi and
            %adapts them for the modified (i.e. reduced) free set. 
        end
        
    end
end