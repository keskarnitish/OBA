function HV = reduced_invHV_lbfgs(V,S,Y,m,k,m_max,F)
HV = zeros(size(V));
if(m==0)
    theta = 1;
    HV(F) = (1/theta)*V(F);
else
    if (k <= m_max)  %Have not filled in memory
        S = S(:,1:m);  Y = Y(:,1:m);
        y = Y(:,m); SY = S'*Y; SS = S'*S;
        theta = (y'*y)/SY(m,m); %initial Hessian approximation is theta*eye(n)
        PW = zeros(length(F),2*m);
        PW(:,1:m)=Y(F,:);  PW(:,m+1:2*m)=theta*S(F,:);
        size_invN = m;
    else %Need to do wrap around from oldest to newest
        %S = [S(:,m+1:m_max), S(:,1:m)];
        %Y = [Y(:,m+1:m_max), Y(:,1:m)];
        y = Y(:,m);  s = S(:,m);
        theta = (y'*y)/(s'*y); %initial Hessian approximation is theta*eye(n)
        S1 = S(:,m+1:m_max); S2 = S(:,1:m);
        Y1 = Y(:,m+1:m_max); Y2 = Y(:,1:m);
        SY = [S1'*Y1, S1'*Y2;
            S2'*Y1, S2'*Y2];
        SS = [S1'*S1, S1'*S2;
            S2'*S1, S2'*S2];
        PW = zeros(length(F),2*m_max);
        PW(:,1:m_max-m)=Y1(F,:);  PW(:,m_max-m+1:m_max)=Y2(F,:);
        PW(:,m_max+1:2*m_max-m)=theta*S1(F,:);  PW(:,2*m_max-m+1:2*m_max)=theta*S2(F,:);
        size_invN = m_max;
    end
    
    %y = Y(:,m);  SY = S'*Y; SS = S'*S;
    D=diag(diag(SY)); L=tril(SY)-D;
    %theta = (y'*y)/SY(m,m); %initial Hessian approximation is theta*eye(n)
    invM = ([-D , L';
        L , theta*SS]);  % M : 2mx2m dimensional
    %PW = [Y(F,:) , theta*S(F,:)];  % PW : |F|x2m dimensional
    invN = (eye(2*size_invN) - (1/theta) * (invM \ (PW' * PW)));  % N : 2m x 2m
    HV(F) = (1/theta)*V(F) + (1/theta)^2 * (PW * (invN \ (invM \ (PW' * V(F)))));
end


end