function V_now=functionv(x,t)
v=zeros(size(x));
V_now1=diag(v,0);
V_now=sparse(V_now1);
end