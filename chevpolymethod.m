function psi_chevpoly=chevpolymethod(psi0,dt,x,J,time_step,dx,delta)
T1=ones(J+1,1);
psi_chevpoly=zeros(J+1,time_step+1);
psi_chevpoly(:,1)=psi0.';
T=spdiags([-T1,2*T1,-T1],-1:1,J+1,J+1)/(2*dx^2);
t_now=0;
for i=1:time_step
    
   
    V=functionv(x,t_now);
    t_now=t_now+dt;
    Hmax=max(max(V))+1/dx^2;
    tau=Hmax*dt;
    H=(T+V)/Hmax;

    Jn_tau=[besselj(0,tau),-2i*besselj(1,tau)];
    ii=3;
    while 1
        tmp=2*(-1i)^(ii-1)*besselj(ii-1,tau);
        if abs(tmp)<delta
            break
        end
        Jn_tau=[Jn_tau,tmp];
        ii=ii+1;
    end
    Npoly=length(Jn_tau);
    psi_chevexp(:,1)=psi_chevpoly(:,i);
    psi_chevexp(:,2)=H*psi_chevexp(:,1);
    for jj=3:Npoly
        psi_chevexp(:,jj)=2*H*psi_chevexp(:,jj-1);
        psi_chevexp(:,jj)=psi_chevexp(:,jj)-psi_chevexp(:,jj-2);
    end
    psi_chevpoly(:,i+1)=psi_chevexp*Jn_tau.';
    

end