function pshi45=myode45(t,dt,psi0,x,J,time_step,dx)
%t_delta=1e-3;
t_delta_change=dt;
iterate_max=1000;
t0=0:t_delta:((iterate_max-1)*t_delta);
psi45=zeros(J+1,time_step+1);
psi45(:,1)=psi0.';
T=spdiags([-T1,2*T1,-T1],-1:1,J+1,J+1)/(2*dx^2);
t_now=0;
%f_n=zeros(12,iterate_max);

s1=myode(0,c_tmp);
t_tol=0;
tol=1e-3;
tmin=1e-6;
y1=[];
for i=1:time_step
     
c_tmp=psi45(:,i);

for iterate=1:1:(iterate_max)

    V1=functionv(x,t_now+t_delta_change/5);
    H1=T+V1;
    V2=functionv(x,t_now+3*t_delta_change/10);
    H2=T+V2;
    V3=functionv(x,t_now+4*t_delta_change/5);
    H3=T+V3;
    V4=functionv(x,t_now+8*t_delta_change/9);
    H4=T+V4;
    V5=functionv(x,t_now+t_delta_change);
    H5=T+V5;
    
        ws1=c_tmp+t_delta_change*s1/5;
        %s2=myode(t_tol+t_delta_change/5,c_tmp+t_delta_change*s1/5);
        s2=-1i*H1*(c_tmp+t_delta_change*s1/5);
        ws2=c_tmp+t_delta_change*(3*s1/40+9*s2/40);
        s3=-1i*H2*(c_tmp+t_delta_change*(3*s1/40+9*s2/40));
        ws3=c_tmp+t_delta_change*(44*s1/45-56*s2/15+32*s3/9);
        s4=-1i*H3*(c_tmp+t_delta_change*(44*s1/45-56*s2/15+32*s3/9));
        ws4=c_tmp+t_delta_change*(19372*s1/6561-25360*s2/2187+64448*s3/6561-212*s4/729);
        s5=-1i*H4*(c_tmp+t_delta_change*(19372*s1/6561-25360*s2/2187+64448*s3/6561-212*s4/729));
        ws5=c_tmp+t_delta_change*(9017*s1/3168-355*s2/33+46732*s3/5247+49*s4/176-5103*s5/18656);
        s6=-1i*H5*(c_tmp+t_delta_change*(9017*s1/3168-355*s2/33+46732*s3/5247+49*s4/176-5103*s5/18656));
        z=c_tmp+t_delta_change*(35*s1/384+500*s3/1113+125*s4/192-2187*s5/6784+11*s6/84);
        s7=-1i*H5*(z);
        normz=norm(z);
    total_err=t_delta_change*norm(71*s1/57600-71*s3/16695+71*s4/1920-17253*s5/339200+s6*22/525-s7/40);
    relative_err=total_err/normz;
    if t_delta_change~=tmin
        if relative_err>tol
            t_delta_change=t_delta_change*0.8*(tol/relative_err)^(1/5);
            if t_delta_change<tmin 
                t_delta_change=tmin;
                iterate=iterate-1;
                continue;
            else
                iterate=iterate-1;
                continue;
            end
        end
    end
    
    
    if (t_tol+t_delta_change)==t_delta
        y1=[y1,z];
        break;
    end
    c_tmp=z;
    
    if (t_tol+t_delta_change)>t_delta
        t_delta_change=t_delta-t_tol;
    end
    t_tol=t_tol+t_delta_change;
    
end
psi45(:,i+1)=z;
t_now=t_now+dt;
end
end
