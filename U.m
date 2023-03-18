function U = U(dt)
    H1 = H(dt/2);
    H2 = H(dt);
    H3 = H(dt/2);
    U = expm(-1i*H2*dt)*expm(-1i*H1*dt)*expm(-1i*H3*dt);
end