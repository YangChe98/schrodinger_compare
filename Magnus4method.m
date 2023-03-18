function psi_ma4=Magnus4method(psi0,dt0,x,J,time_step,dx)
T1=ones(J+1,1);
psi_ma4=zeros(J+1,time_step+1);
psi_ma4(:,1)=psi0.';
T=spdiags([-T1,2*T1,-T1],-1:1,J+1,J+1)/(2*dx^2);
t_now=0;
% v=zeros(size(x));
% V_now1=diag(v,0);
% V_now=sparse(V_now1);
% V_now=functionv(x,t_now);
sparsei = speye(J+1);
for i=1:time_step
    t_forward=t_now;
    c1=1/2-sqrt(3)/6;
    c2=1/2+sqrt(3)/6;
    V1=functionv(x,t_now+dt0*c1);
    V2=functionv(x,t_now+dt0*c2);
    A1=-1i*(T+V1);
    A2=-1i*(T+V2);
    Omega=dt0*(A1+A2)/2-sqrt(3)*dt0^2*(A1*A2-A2*A1)/12;
    t_now=t_now+dt0;
    
    
    %psi_now=U2\(U1*psi_cn(:,i));
    psi_ma4(:,i+1)=expm(Omega)*psi_ma4(:,i);
end
end