% Section 4: The Impact Map

 function Delta=impact_map(x,data)
    % extract variables and functions from input arguments
    q = x(1:end/2);
    qdot = x(end/2+1:end);
    Dbarfun = data.Dbar;
    Efun = data.E;

    delta1_qdot = [eye(5) zeros(5,4)]*( [ Dbarfun(q) -Efun(q).' ; Efun(q) zeros(2,2) ]\[Dbarfun(q)*[eye(5); zeros(2,5)]; zeros(2,5)] ) * qdot;

    delta2_R = [1  1  1  1 0;...
                0  0  0 -1 0; ...
                0  0 -1  0 0; ...
                0 -1  0  0 0; ...
                0  0 -1  0 1];
    delta2_d = [-3; 2; 2; 2; 1].*pi;

    q_final = delta2_R*q + delta2_d;
    qdot_final = delta2_R*delta1_qdot;
    Delta = [q_final; qdot_final];
 end