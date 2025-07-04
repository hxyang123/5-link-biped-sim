% Create an ODE function biped.m accepting the structure array data. The function declaration
% should be this: function xdot=biped(t,x,data). Here, x is the robot state, [q;  ̇q] and xdot is its
% time derivative. The function therefore implements the closed-loop vector field.
% Within this function, extract from x the subvectors q and qdot, and extract from the structure
% data the variables you need for the ODE. You need to compute τ as in (2) and compute  ̈q =
% D−1(q)(−C(q,  ̇q)  ̇q − ∇q P(q) + Bτ). Once you’ve done that, you’ll set xdot = [qdot;qddot].

function xdot=biped(t,x,data)
    % extract variables and functions needed for the ODE
    q = x(1:end/2);
    qdot = x(end/2+1:end);
    Dfun = data.D;
    Cfun = data.C;
    Gfun = data.G;
    B = data.B;
    Kp = data.Kp;
    Kd = data.Kd;

    phi2 = data.phi2;
    phiprime2 = data.phiprime2;
    phipprime2 = data.phipprime2;
    
    phi3 = data.phi3;
    phiprime3 = data.phiprime3;
    phipprime3 = data.phipprime3;

    phi4 = data.phi4;
    phiprime4 = data.phiprime4;
    phipprime4 = data.phipprime4;

    phi5 = data.phi5;
    phiprime5 = data.phiprime5;
    phipprime5 = data.phipprime5;

    % added info
    q_tilde = data.q_tilde;
    q_plus = data.q_plus;
    q_minus = data.q_minus;

    H = [phiprime2((q_plus(1) - q(1))/q_tilde(1))/q_tilde(1) 1 0 0 0;
        phiprime3((q_plus(1) - q(1))/q_tilde(1))/q_tilde(1) 0 1 0 0;
        phiprime4((q_plus(1) - q(1))/q_tilde(1))/q_tilde(1) 0 0 1 0;
        phiprime5((q_plus(1) - q(1))/q_tilde(1))/q_tilde(1) 0 0 0 1];

    y = [q(2) - phi2((q_plus(1) - q(1))/q_tilde(1));
        q(3) - phi3((q_plus(1) - q(1))/q_tilde(1));
        q(4) - phi4((q_plus(1) - q(1))/q_tilde(1));
        q(5) - phi5((q_plus(1) - q(1))/q_tilde(1))];
    ydot = [qdot(2) + phiprime2((q_plus(1) - q(1))/q_tilde(1))*qdot(1)/q_tilde(1);
        qdot(3) + phiprime3((q_plus(1) - q(1))/q_tilde(1))*qdot(1)/q_tilde(1);
        qdot(4) + phiprime4((q_plus(1) - q(1))/q_tilde(1))*qdot(1)/q_tilde(1);
        qdot(5) + phiprime5((q_plus(1) - q(1))/q_tilde(1))*qdot(1)/q_tilde(1)];

    H2 = [-phipprime2((q_plus(1) - q(1))/q_tilde(1))/q_tilde(1)^2;
          -phipprime3((q_plus(1) - q(1))/q_tilde(1))/q_tilde(1)^2;
          -phipprime4((q_plus(1) - q(1))/q_tilde(1))/q_tilde(1)^2;
          -phipprime5((q_plus(1) - q(1))/q_tilde(1))/q_tilde(1)^2];


    % matlab warning: INV(A)*b can be slower and less accurate than A\b.
    % Consider using A\b for INV(A)*b
    % tau = inv( H*inv( Dfun(q) )*B )*( H*inv( Dfun(q) )*( Cfun(q,qdot)*qdot + Gfun(q) ) - Kp*sin.(H*q-qref) - Kd*H*qdot  );
    tau = ( H*( Dfun(q)\B ) )\( H*( Dfun(q) \ ( Cfun(q,qdot)*qdot + Gfun(q) ) ) -H2- Kp*sin(y) - Kd*ydot  );
    % qddot = inv( Dfun(q) )*( -Cfun(q,qdot)*qdot - Gfun(q) + B*tau );
    qddot = Dfun(q) \ ( -Cfun(q,qdot)*qdot - Gfun(q) + B*tau );

    xdot= [qdot; qddot];

end
