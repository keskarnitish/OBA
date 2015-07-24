function [PW,theta,invM] = lbfgs(S,Y,m,k,m_max)
% B = \theta I - WMW'
%Page 4 of Byrd, La, Nocedal, Zha
if(m==0)
    theta = 1;
    PW = [];  invM = [];
else
    %         if (k <= m_max)  %Have not filled in memory
    S = S(:,1:m);  Y = Y(:,1:m);
    y = Y(:,m); SY = S'*Y; SS = S'*S;
    theta = (y'*y)/SY(m,m); %initial Hessian approximation is theta*eye(n)
    PW = zeros(length(y),2*m);
    PW(:,1:m)=theta*S;  PW(:,m+1:2*m)=Y;
    
    %         else %Need to do wrap around from oldest to newest
    %             %S = [S(:,m+1:m_max), S(:,1:m)];
    %             %Y = [Y(:,m+1:m_max), Y(:,1:m)];
    %             y = Y(:,m);  s = S(:,m);
    %             theta = (y'*y)/(s'*y); %initial Hessian approximation is theta*eye(n)
    %             S1 = S(:,m+1:m_max); S2 = S(:,1:m);
    %             Y1 = Y(:,m+1:m_max); Y2 = Y(:,1:m);
    %             SY = [S1'*Y1, S1'*Y2;
    %                 S2'*Y1, S2'*Y2];
    %             SS = [S1'*S1, S1'*S2;
    %                 S2'*S1, S2'*S2];
    %             %theta = (y'*y)/SY(m,m); %initial Hessian approximation is theta*eye(n)
    %             PW = zeros(length(y),2*m_max);
    %             PW(:,1:m_max-m)=theta*S1;  PW(:,m_max-m+1:m_max)=theta*S2;
    %             PW(:,m_max+1:2*m_max-m)=Y1;  PW(:,2*m_max-m+1:2*m_max)=Y2;
    %         end
    
    %y = Y(:,m);  SY = S'*Y; SS = S'*S;
    D=diag(diag(SY)); L=tril(SY)-D;
    %theta = (y'*y)/SY(m,m); %initial Hessian approximation is theta*eye(n)
    invM = ([theta*SS , L;
        L' , -D]);  % M : 2mx2m dimensional
    %PW = [theta*S, Y];  % PW : nx2m dimensional
    %PWV = PW'*V;
    %HV = theta*V - PW*(invM\PWV);
end

end
