clear;
close all;
%% 1) Define numerical variables

% link lengths (m)
l1 = 0.5;
l2 = 0.5;
l3 = 0.3;

%masses (Kg)
m1 = 0.05;
m2 = 0.5;
m3 = 0.3;
m4 = 0.5;
m5 = 0.05;
m6 = 0.5;

% gravity
g = 9.81; %m/s^2

% degrees of freedom
n = 5;

% 2) Define real symbolic variables
syms t q1 q2 q3 q4 q5 q1dot q2dot q3dot q4dot q5dot real
syms x1 x2 x1dot x2dot real
syms tau1 tau2 tau3 tau4 real

q = [q1; q2; q3; q4; q5];
qdot = [q1dot; q2dot; q3dot; q4dot; q5dot];
tau = [tau1; tau2; tau3; tau4];
x = [x1; x2];
xdot = [x1dot; x2dot];
qbar = [q; x];
qbardot = [qdot; xdot];

% 3) Find total kinetic energy for unpinned robot
%%% a) positions of six masses [ri = ....]
% Note - there should also be a handwritten submission for this derivation
% Note - assumes stance foot is at x 

% m1
r1 = x;
% r2
r2 = r1 + l1 * [cos(q1); sin(q1)];
% r3
r3 = r2 + l2 * [cos(q1+q2); sin(q1+q2)];
% r6
r6 = r3 + l3 * [cos(q1+q2+q5); sin(q1+q2+q5)];

% r4
r4 = r3 + l2 * [cos(q1+q2+q3); sin(q1+q2+q3)];
% r5
r5 = r4 + l1 * [cos(q1+q2+q3+q4); sin(q1+q2+q3+q4)];

%%% b)  Symbolic expressions for r and rdot

r1dot = jacobian(r1, qbar)*qbardot;
r2dot = jacobian(r2, qbar)*qbardot;
r3dot = jacobian(r3, qbar)*qbardot;
r4dot = jacobian(r4, qbar)*qbardot;
r5dot = jacobian(r5, qbar)*qbardot;
r6dot = jacobian(r6, qbar)*qbardot;
%

%%% c) Kinetic energy
K = 0.5 * (m1 * (r1dot'*r1dot) + m2 * (r2dot'*r2dot) ...
    + m3 * (r3dot'*r3dot) + m4 * (r4dot'*r4dot) ...
    + m5 * (r5dot'*r5dot) + m6 * (r6dot'*r6dot));
K = simplify(K);
%
%%% d) Dbar
% Dbar=simplify(hessian(K,qbardot));
Dbar = hessian(K,qbardot);


% 
%%% e) D matrix 
% This matrix D(q) is the kinetic energy of the pinned robot
% Submatrix of Dbar; first 5 rows and columns

D = Dbar(1:5, 1:5); 

% 4
% Write the total potential energy of the pinned robot in terms of q and create a variable P containing
% it. Create a column vector G containing the symbolic gradient (i.e., the transposed Jacobian) of P,
h1 = 0;
h2 = l1*sin(q1);
h3 = h2 + l2*sin(q1+q2);
h4 = h3 + l2*sin(q1+q2+q3);
h5 = h4 + l1*sin(q1+q2+q3+q4);
h6 = h3 + l3*sin(q1+q2+q5);

h_mat = [h1,h2, h3, h4, h5, h6];
M_mat = [m1, m2, m3, m4, m5, m6];

P = g*dot(M_mat, h_mat);
G = jacobian(P,q)';

% 5
% Define the input matrix B of the robot, of dimension 5 × 4. Note that the first configuration variable,
% q1 in unactuated, while q2, . . . , q5 are directly actuated.

B = [zeros(1, 4); eye(4,4);];
%6
% Using symbolic differentiation and the formula given in class for the Christoffel coefficients of D
% (note: D, not Dbar), find the Coriolis matrix C(q, ˙q) and name it C.

C = sym(zeros(n,n));
for i=1:n
    for j=1:n
        for k=1:n
            C(k,j)=C(k,j)+1/2*(diff(D(k,j),q(i))+diff(D(k,i),q(j))-diff(D(i,j),q(k)))*qdot(i);
        end
    end
end

% 7 
% The impact map ∆ requires Dbar(q), which you’ve already computed, and the matrix-valued function
% E(q) = [(dpx)q I_2]
% where (dpx)q is the Jacobian of px(q), the vector with head at the swing foot and tail at x (see
% Figure 1). Define symbolic variables px and E containing px(q) and the expression above for E(q).

px = r5 - r1;
E=[jacobian(px,q), eye(2)];


% 8
% Turn the symbolic expressions in D,Dbar,px,E,C,G into Matlab functions with appropriate inputs.

Dfun=matlabFunction(D,'Vars',{q});
Dbarfun = matlabFunction(Dbar, 'Vars', {qbar});
pxfun = matlabFunction(px, 'Vars', {q});
Efun = matlabFunction(E, 'Vars', {q});
Cfun = matlabFunction(C, 'Vars',{q, qdot});
Gfun = matlabFunction(G, 'Vars',{q});

% 9)
% Create a structure array named data containing these objects: Dfun,Dbarfun,Efun,Cfun,Gfun,B.
% For instance, the command data.D=Dfun will create a field D in the structure containing the function
% Dfun. Later on, you will pass this structure array to the ode45 function for numerical integration,
% and to the impact_map function for computation of the impact map.

data = struct;
data.D = Dfun;
data.Dbar = Dbarfun;
data.E = Efun;
data.C = Cfun;
data.G = Gfun;
data.B = B;

% 1) Define numerical value

% Define two control gains Kp,Kd for the controller (2)
% placing the roots of the polynomial λ2 + Kdλ + Kp at {−20, −20}
Kd = 40;
Kp = 400;
H = [zeros(4,1) eye(4)];

data.H = H;
data.Kp = Kp;
data.Kd = Kd;

%% [NEW] Adding the Creative section
% The outcome of this part should be a parameter vector a, whose components
% a_1,ldots,a_k define the polynomial phi_a(theta)=a_1 + a_2 theta + ... +
% a_k theta^(k-1)

fprintf('\n Determining VHC...\n')

% select leg angle of aperture, extract matrix I
alpha = pi - pi/6;
beta = 0.316637; % (0, pi)
v1 = -0.894373;
v2 = 1.9;

% Define numerical qref
q2ref = pi- alpha;
q4ref = pi+alpha; 
q5ref = -pi/4;

% set up equation

% set up q_plus and q_minus and verify with impact map
q_minus = [(pi-beta)/2 - (pi-alpha)/2; 
            q2ref;
            pi + beta;
            q4ref;
            q5ref+beta];
q_plus = [(beta+pi)/2 - (pi-alpha)/2; 
            q2ref;
            pi - beta;
            q4ref;
            q5ref;];
q_tilde = q_plus - q_minus;

% from acrobot simple VHC
% I_delta_num = [-0.6864    0.3621;
%    -0.1459   -0.1081];
% f_v1 = 0.4246;

% Calculating I and f(v1) from impact_map
q_input = q_minus;
delta1_qdot = [eye(5) zeros(5,4)]*( [ Dbarfun(q_input) -Efun(q_input).' ; Efun(q_input) zeros(2,2) ]\...
                [Dbarfun(q_input)*[eye(5); zeros(2,5)]; zeros(2,5)] ) * qdot;
delta1_qdot = simplify(delta1_qdot);
delta2_R = [1  1  1  1 0;...
            0  0  0 -1 0; ...
            0  0 -1  0 0; ...
            0 -1  0  0 0; ...
            0  0 -1  0 1];
qdot_final = delta2_R*delta1_qdot;
I_delta = double(jacobian(qdot_final,qdot));
f_v1 = -q_tilde(1)*(-I_delta(3,1)*q_tilde(1) + I_delta(3,3)*v1)/...
        (-I_delta(1,1)*q_tilde(1) + I_delta(1,3)*v1);

% VHC for q2
A_VHC2 = [1 0 0 0 0 0;
         1 1 1 1 1 1;
         1 0.5 0.5^2 0.5^3 0.5^4 0.5^5];
b_VHC2 = [q_plus(2); q_minus(2);q2ref];

% VHC for q3
A_VHC3 = [1 0 0 0 0 0;
         1 1 1 1 1 1;
         1 0.5 0.5^2 0.5^3 0.5^4 0.5^5;
         0 1 0 0 0 0;
         0 1 2 3 4 5;
         0 1 2*0.5 3*0.5^2 4*0.5^3 5*0.5^4;];
b_VHC3 = [q_plus(3); q_minus(3); pi; f_v1; v1; v2];

% VHC for q4
A_VHC4 = [1 0 0 0 0 0;
         1 1 1 1 1 1;
         1 0.5 0.5^2 0.5^3 0.5^4 0.5^5];
b_VHC4 = [q_plus(4); q_minus(4); q4ref-pi/4];

% VHC for q5
A_VHC5 = [1 0 0 0 0 0;
         1 1 1 1 1 1;
         1 0.5 0.5^2 0.5^3 0.5^4 0.5^5];
b_VHC5 = [q_plus(5); q_minus(5); q5ref+beta/2];



% solve for a
a2 = linsolve(A_VHC2,b_VHC2);
a3 = linsolve(A_VHC3,b_VHC3);
a4 = linsolve(A_VHC4,b_VHC4);
a5 = linsolve(A_VHC5,b_VHC5);

% HERE WE DEFINE THE FUNCTION phi_a AND ITS DERIVATIVES
fprintf('\n Define phi_a and its derivatives...\n')
a2=flip(a2);
a3=flip(a3);
a4=flip(a4);
a5=flip(a5);

% define all phi and dericatives
phi2=@(theta) polyval(a2,theta);
phiprime2=@(theta) polyval(polyder(a2),theta);
phipprime2=@(theta) polyval(polyder(polyder(a2)),theta);

phi3=@(theta) polyval(a3,theta);
phiprime3=@(theta) polyval(polyder(a3),theta);
phipprime3=@(theta) polyval(polyder(polyder(a3)),theta);

phi4=@(theta) polyval(a4,theta);
phiprime4=@(theta) polyval(polyder(a4),theta);
phipprime4=@(theta) polyval(polyder(polyder(a4)),theta);

phi5=@(theta) polyval(a5,theta);
phiprime5=@(theta) polyval(polyder(a5),theta);
phipprime5=@(theta) polyval(polyder(polyder(a5)),theta);

% Using phi and its derivatives, below you should define functions sigma,
% sigmaprime, sigmapprime.

sigma = @(theta) sigma_fun(theta,q_plus,q_tilde,phi2,phi3,phi4,phi5,alpha);
sigmaprime = @(theta) sigmaprime_fun(theta,q_tilde,phiprime2,phiprime3,phiprime4,phiprime5);
sigmapprime = @(theta) sigmapprime_fun(theta,phipprime2,phipprime3,phipprime4,phipprime5);

% Add information for the datastructure
data.phi2=phi2;
data.phiprime2=phiprime2;
data.phipprime2=phipprime2;

data.phi3=phi3;
data.phiprime3=phiprime3;
data.phipprime3=phipprime3;

data.phi4=phi4;
data.phiprime4=phiprime4;
data.phipprime4=phipprime4;

data.phi5=phi5;
data.phiprime5=phiprime5;
data.phipprime5=phipprime5;

data.sigma=sigma;
data.sigmaprime=sigmaprime;
data.sigmapprime=sigmapprime;

% added info
data.q_tilde = q_tilde;
data.q_plus = q_plus;
data.q_minus = q_minus;

% HERE WRITE CODE TO TEST WHETHER YOUR VHC WORKS
fprintf('\n Validating VHC...\n')

theta = 0:0.001:1;
q_test = [];
% Plot q and check if it passes q plus q minus and q bar
for i = 1:1001
    q_test = [q_test sigma(theta(i))];
end

figure(2)
hold on

yyaxis right
plot(q_test(1,:), q_test(2,:),'--');
plot(q_test(1,:), q_test(4,:),':');
plot(q_test(1,:), q_test(5,:),'-.');
yyaxis left
plot(q_test(1,:), q_test(3,:));
xlabel 'q1'
ylabel 'q3'

% check is q3 passes pi when p1 passes alpha/2
title 'VHC q = sigma(theta)'
plot(q_test(1,1), q_test(3,1),"o");
plot(q_test(1,end-1), q_test(3,end-1),"o");
plot(alpha/2,pi,"o");
text(q_test(1,1), q_test(3,1),"   q plus");
text(q_test(1,end-1), q_test(3,end-1),"   q minus");
text(alpha/2,pi,"   q bar");

% verify that the curve is a regular VHC
B_perp = null(B');
regularVHC = [];
for i = 1:1001
    regularVHC = [regularVHC; B_perp.'*Dfun(q_test(:,i))*sigmaprime(theta(i))];
end
figure(3)
plot(theta,regularVHC);
title 'Regular VHC check'

% verify whether the constraints induces a stable limit cycle
% Building psi
Psi1 = [];
Psi2 = [];
for i = 1:1001
    denominator = B_perp.'*Dfun(q_test(:,i))*sigmaprime(theta(i));
    Psi1= [Psi1; -B_perp.'*Gfun(q_test(:,i))/denominator];
    Psi2= [Psi2; -B_perp.'*(Dfun(q_test(:,i))*sigmapprime(theta(i))+Cfun(q_test(:,i),sigmaprime(theta(i)))*sigmaprime(theta(i)))/denominator];
end

% calculate mass, velocity, and delta
M_theta_num = exp(-2*cumtrapz(theta,Psi2));
V_theta_num = -cumtrapz(theta,Psi1.*M_theta_num);

M_theta = @(theta) M_theta_num((theta*1000)+1);
V_theta = @(theta) V_theta_num((theta*1000)+1);
q_input = sigma(1);
delta1_qdot = [eye(5) zeros(5,4)]*( [ Dbarfun(q_input) -Efun(q_input).' ; Efun(q_input) zeros(2,2) ]\...
                [Dbarfun(q_input)*[eye(5); zeros(2,5)]; zeros(2,5)] ) * qdot;
delta1_qdot = simplify(delta1_qdot);

delta2_R = [1  1  1  1 0;...
            0  0  0 -1 0; ...
            0  0 -1  0 0; ...
            0 -1  0  0 0; ...
            0  0 -1  0 1];

qdot_final = delta2_R*delta1_qdot;
I_delta = double(jacobian(qdot_final,qdot));

delta_theta = dot(sigmaprime(0), I_delta*sigmaprime(1))/...
              (sigmaprime(0).'*sigmaprime(0));
% stable hybrid limit cycle check
if delta_theta^2/M_theta(1) > 0 && delta_theta^2/M_theta(1) < 1
    fprintf(" Mass Check Passed\n")
else
    fprintf(" Mass Check Failed\n")
end

if V_theta(1)*delta_theta^2/(M_theta(1)-delta_theta^2) + max(V_theta_num) < 0
    fprintf(" Velocity Check Passed\n")
else
    fprintf(" Velocity Check Failed\n")
end
 %% Section 5: Numerical Simulation

ops= odeset('reltol',1e-7,'abstol',1e-7,'Events',@ground_impact);
dt=1/60; % 60 fps; time increment in simulations and animations

fprintf('\n Simulating...\n')

simtype = 1;

if simtype == 1 % On the limit cycle
    % from acrobot VHC
    % q0 = 1.7291 2.8250 qdot0 = -2.1172 2.8393
    q0 = sigma(0);
    qdot0 = [-2.1172;
            0;
            2.8393;
            0;
            0];

elseif simtype == 2 % On the constraint manifold
    q0 = sigma(0.7);
    qdot0 = sigmaprime(0.7)*delta_theta*sqrt(-2*V_theta(1)/(M_theta(1)-delta_theta^2));

else % outside of constraint manifold
    q0_eps = 0.01;
    qdot0_eps = 0.02;
    q0 = sigma(0.7) + q0_eps;
    qdot0 = sigmaprime(0.7)*delta_theta*sqrt(-2*V_theta(1)/(M_theta(1)-delta_theta^2)) + qdot0_eps;
end

% initial condition

T=[];
X=[];
Te=[];
Ie=[];
Xe=[];
post_impact_state=[q0;qdot0];
% Simulate number_steps steps
for step=1:25
    fprintf('\n...step %d...\n',step);
    [t,x,te,xe,ie]=ode45(@(t,x) biped(t,x,data),0:dt:10,post_impact_state,ops);
    % Application of the impact map
    impact_state=xe(end,:)';
    post_impact_state=impact_map(impact_state,data);
    T{step}=t;
    X{step}=x;
    Ie{step}=ie;
    Te{step}=te;
    Xe{step}=xe;
end

fprintf('\n Setting up animation...\n')



%% Animation of the simulation results
figure(1);
ref=0;time_passed=0;step=1;
Axis=[-1 4 0 2];
Time=text(-1+2,1.8,['time= ','0',' secs,',' step= ',num2str(step)]);
axis(Axis);
q = q0;
q1=q(1);q2=q(2);q3=q(3);q4=q(4);q5=q(5);
xdata=0;
ydata=0;
l=[l1 l2 l2 l1 l3];
Q=[q1 q1+q2 q1+q2+q3 q1+q2+q3+q4 q1+q2+q5];
for j=1:4
    xdata=[xdata xdata(end)+l(j)*cos(Q(j))];
    ydata=[ydata ydata(end)+l(j)*sin(Q(j))];
end
xdata=[xdata xdata(3)+l(5)*cos(Q(5))];
ydata=[ydata ydata(3)+l(5)*sin(Q(5))];

link1=line([xdata(1) xdata(2)],[ydata(1) ydata(2)],'color','red','linewidth',2);
link2=line([xdata(2) xdata(3)],[ydata(2) ydata(3)],'color','red','linewidth',2);
link3=line([xdata(3) xdata(4)],[ydata(3) ydata(4)],'linewidth',2);
link4=line([xdata(4) xdata(5)],[ydata(4) ydata(5)],'linewidth',2);
link5=line([xdata(3) xdata(6)],[ydata(3) ydata(6)],'linewidth',2);

fprintf('\n Animation is ready...\n')
ref=0; % This variable keeps track of the position of the stance foot accross multiple steps

% v = VideoWriter("P3_5link_robot");
% open(v)

animation_slowdown_factor=2; % >1 means slow down
for step=1:length(Ie)
    t=T{step};
    x=X{step};
    xe=Xe{step};
    xe=xe(end,:);
    x_data_5 = 0;
    for k=2:length(t)
        t0=clock;
        drawnow;
        q=x(k,1:5)';
        q1=q(1);q2=q(2);q3=q(3);q4=q(4);q5=q(5);
        Q=[q1 q1+q2 q1+q2+q3 q1+q2+q3+q4 q1+q2+q5];
        xdata=0;
        ydata=0;
        for j=1:4
            xdata=[xdata xdata(end)+l(j)*cos(Q(j))];
            ydata=[ydata ydata(end)+l(j)*sin(Q(j))];
        end
        xdata=ref+[xdata xdata(3)+l(5)*cos(Q(5))];
        ydata=[ydata ydata(3)+l(5)*sin(Q(5))];
        
        x_data_5 = xdata(5);

        set(link1,'xdata',[xdata(1) xdata(2)],'ydata',[ydata(1) ydata(2)]);
        set(link2,'xdata',[xdata(2) xdata(3)],'ydata',[ydata(2) ydata(3)]);
        set(link3,'xdata',[xdata(3) xdata(4)],'ydata',[ydata(3) ydata(4)]);
        set(link4,'xdata',[xdata(4) xdata(5)],'ydata',[ydata(4) ydata(5)]);
        set(link5,'xdata',[xdata(3) xdata(6)],'ydata',[ydata(3) ydata(6)]);
        set(Time,'String',['time= ',num2str(round(time_passed+t(k),1)),' secs,',' step= ',num2str(step)]);
        current_axis=gca;
        if ref>.95*current_axis.XLim(end)
            current_axis.XLim=[.95*ref .95*ref+5];
            Time.Position=[.95*ref+2 1.8 0];
            Axis=axis;
        else
            axis(Axis)
        end
        while etime(clock,t0)<animation_slowdown_factor*(t(k)-t(k-1))
        end
        % frame = getframe(gcf);
        % writeVideo(v,frame)
    end
    time_passed=time_passed+t(end);
    % ref=ref+l*(cos(xe(1))+cos(xe(1)+xe(2)));
    ref=x_data_5;
end
% close(v);
%%

% Define sigma functions
function sigma = sigma_fun(theta,q_plus,q_tilde,phi2,phi3,phi4,phi5,alpha)
    sigma = q_plus;
    sigma(1) = q_plus(1) - theta*q_tilde(1);
    sigma(2) = phi2(theta);
    sigma(3) = phi3(theta);
    sigma(4) = phi4(theta);
    sigma(5) = phi5(theta);
end

function sigmaprime = sigmaprime_fun(theta,q_tilde,phiprime2,phiprime3,phiprime4,phiprime5)
    sigmaprime = zeros(5,1);
    sigmaprime(1) = -q_tilde(1);
    sigmaprime(2) = phiprime2(theta);
    sigmaprime(3) = phiprime3(theta);
    sigmaprime(4) = phiprime4(theta);
    sigmaprime(5) = phiprime5(theta);
    
end

function sigmapprime = sigmapprime_fun(theta,phipprime2,phipprime3,phipprime4,phipprime5)
    sigmapprime = zeros(5,1);
    sigmapprime(2) = phipprime2(theta);
    sigmapprime(3) = phipprime3(theta);
    sigmapprime(4) = phipprime4(theta);
    sigmapprime(5) = phipprime5(theta);
end