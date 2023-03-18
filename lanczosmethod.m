function psi_lan=lanczosmethod(psi0,dt,x,J,time_step,dx,M)
T1=ones(J+1,1);
psi_lan=zeros(J+1,time_step+1);
psi_lan(:,1)=psi0.';
T=spdiags([-T1,2*T1,-T1],-1:1,J+1,J+1)/(2*dx^2);
t_now=0;
V_lan=zeros(J+1,M);
H_lan=zeros(M,M);
e1=zeros(M,1);
e1(1,1)=1;
for i=1:time_step
    V=functionv(x,t_now);
     t_now=t_now+dt;
    H=T+V;
    V_lan(:,1)=psi_lan(:,i)/sqrt(trapz(x,abs(psi_lan(:,i)).^2));
    H_lan(1,1)=trapz(x, conj(V_lan(:,1)).*(H*V_lan(:,1)));

    for k=1:(M-1)
        V_lan(:,k+1)=H*(V_lan(:,k));
        for j=1:(k+1)
            if k-j<=1
               V_lan(:,k+1)=V_lan(:,k+1)-H_lan(j,k)*V_lan(:,j);
            end
        end
        V_lan(:,k+1)=V_lan(:,k+1)/sqrt(trapz(x,abs(V_lan(:,k+1)).^2));

        H_lan(k+1,k+1)=trapz(x,conj(V_lan(:,k+1)).*(H*V_lan(:,k+1)));
        H_lan(k+1,k)=trapz(x,conj(V_lan(:,k+1)).*(H*V_lan(:,k)));
        H_lan(k,k+1)=conj(H_lan(k+1,k));
    end
%         [U,E]=eig(H_lan);

       % phi_lan_tmp=U*exp(-1i*dt*E)*U'*e1;
        phi_lan_tmp=expm(-1i*dt*H_lan)*e1;

        psi_lan(:,i+1)=V_lan* phi_lan_tmp;
    end
end
