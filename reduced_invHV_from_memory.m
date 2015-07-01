function HV = reduced_invHV_from_memory(V,PW,theta,invM,F)
HV = zeros(size(V));
if(isempty(PW))
    HV(F) = (1/theta)*V(F);
else
    size_invN2 = size(invM,1); 
    invN = (eye(size_invN2) - (1/theta) * (invM \ (PW' * PW)));  % N : 2m x 2m
    PWF = PW(F,:);
    HV(F) = (1/theta)*V(F) + (1/theta)^2 * (PWF * (invN \ (invM \ (PWF' * V(F)))));
end


end