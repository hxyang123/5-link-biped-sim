%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ECE1658
%% FINAL PROJECT
%% VHC for walking acrobot
%% PARTIAL CODE TO GET STUDENTS STARTED ON THE PROJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

number_steps=25;
symbolic_computations=1;
%% Physical parameters
l=1;
lc=0.5;
m=1;
Iz=1/12*m*l^2;
g=9.81;
incline_degrees = 2;
psi=deg2rad(incline_degrees);
%% Control parameters
Kp=20^2;
Kd=20*2;

%% Symbolic computations
if symbolic_computations
    % Define symbolic variables
    fprintf('\n Initializing symbolic variables...\n')
    % syms l lc Iz m g real
    syms t q1 q2 x1 x2 q1dot q2dot x1dot x2dot tau real
    
    q=[q1;q2];x=[x1;x2]; qbar=[q;x];
    qdot=[q1dot;q2dot];  xdot=[x1dot;x2dot]; qbardot=[qdot;xdot];
    
    fprintf('\n Symbolic model computation...\n')
    
    % Define centres of mass of two links
    rc1=x+lc*[cos(q1);sin(q1)];
    rc2=x+l*[cos(q1);sin(q1)]+lc*[cos(q1+q2);sin(q1+q2)];
    
    % Compute time derivatives of centres of mass
    rc1dot=jacobian(rc1,qbar)*qbardot;
    rc2dot=jacobian(rc2,qbar)*qbardot;
    
    % Define the total kinetic energy of the robot
    K=1/2*m*(rc1dot'*rc1dot+rc2dot'*rc2dot)+1/2*Iz*(q1dot^2+(q1dot+q2dot)^2);
    K=simplify(K);
    
    % Extract the square symmetric matrix of the kinetic energy
    Dbar=simplify(hessian(K,qbardot));
    
    % Extract the matrix of the kinetic energy of the pinned robot
    D = Dbar(1:2,1:2);
    
    % Define the potential energy of the pinnedrobot
    P = m*g*(lc*sin(q1)+l*sin(q1)+lc*sin(q1+q2));
    psi=deg2rad(incline_degrees);
    gvector=[cos(psi) -sin(psi);sin(psi) cos(psi)]*[0; -1];
    P=-(m*g*lc*[cos(q1);sin(q1)]+m*g*[l*cos(q1)+lc*cos(q1+q2);l*sin(q1)+lc*sin(q1+q2)])'*gvector;
    P=simplify(P);
    % Input matrix of pinned robot
    B=[0;1];
    
    % Computation of matrix C(q,qdot) of pinned robot
    C = sym(zeros(2,2));
    for i=1:2
        for j=1:2
            for k=1:2
                C(k,j)=C(k,j)+1/2*(diff(D(k,j),q(i))+diff(D(k,i),q(j))-diff(D(i,j),q(k)))*qdot(i);
            end
        end
    end
    
    % Computation of gradient of the potential for the pinned robot
    
    G = jacobian(P,q)';
    
    % Computation of qddot (pinned robot)
    
    %     qddot = simplify(inv(D)*(-C*qdot-G+B*tau));
    Dfun=matlabFunction(D,'Vars',{q});
    Gfun=matlabFunction(G,'Vars',{q});
    Cfun=matlabFunction(C,'Vars',{[q;qdot]});
    fprintf('\n Impact map computation...\n')
    
    %% Impact map
    
    % Relabelling map
    
    T=-eye(2)*q + [pi;2*pi];
    Delta2=[T;-eye(2)*qdot];
    
    % First impact map
    
    px=l*[cos(q1)+cos(q1+q2);sin(q1)+sin(q1+q2)];
    xo=sym(zeros(2,1));
    
    E=[jacobian(px,q), eye(2)];
    Deltaqdot=[eye(2),zeros(2,4)]*inv([Dbar -E';E zeros(2,2)])*...
        [Dbar*[eye(2);jacobian(xo,q)];zeros(2,2)];
    Delta1=eval([q;Deltaqdot*qdot]);
    Delta1=simplify(Delta1);
    
    % Composition of two maps
    Delta=simplify(subs(Delta2,[q;qdot],Delta1));
    Deltafun=matlabFunction(Delta,'Vars',{[q;qdot]});
    save('walking_acrobot_model','D','Dfun','Cfun','Gfun','Deltafun','B');
else
    fprintf('\nLoading robot model...\n')
    load('walking_acrobot_model');
end
%% HERE WRITE YOUR CODE FOR THE VHC DESIGN
% The outcome of this part should be a parameter vector a, whose components
% a_1,ldots,a_k define the polynomial phi_a(theta)=a_1 + a_2 theta + ... +
% a_k theta^(k-1)

fprintf('\n Determining VHC...\n')

% select leg angle of aperture, extract matrix I
beta = 0.316637; % (0, pi)
v1 = -0.894373;
v2 = 1.9;

% find useful information about the robot model
q_minus = [(pi-beta)/2; pi + beta];
q_plus = [(beta+pi)/2; pi - beta];
Delta_jacobian = jacobian(Delta,[q1dot;q2dot]);
I_delta = Delta_jacobian(3:4,1:2);
I_deltafun=matlabFunction(I_delta,'Vars',{q});
I_delta_num = I_deltafun(q_minus);

% set up equation
k = 6;
q_tilde = q_plus - q_minus;
f_v1 = -q_tilde(1)*(-I_delta_num(2,1)*q_tilde(1) + I_delta_num(2,2)*v1)/...
        (-I_delta_num(1,1)*q_tilde(1) + I_delta_num(1,2)*v1);

% set up equation
A_VHC = [1 0 0 0 0 0;
         1 1 1 1 1 1;
         1 0.5 0.5^2 0.5^3 0.5^4 0.5^5;
         0 1 0 0 0 0;
         0 1 2 3 4 5;
         0 1 2*0.5 3*0.5^2 4*0.5^3 5*0.5^4;];
b_VHC = [q_plus(2); q_minus(2); pi; f_v1; v1; v2];

% solve for a
a = linsolve(A_VHC,b_VHC);

%% HERE WE DEFINE THE FUNCTION phi_a AND ITS DERIVATIVES
fprintf('\n Define phi_a and its derivatives...\n')
a=flip(a);
phi=@(theta) polyval(a,theta);
phiprime=@(theta) polyval(polyder(a),theta);
phipprime=@(theta) polyval(polyder(polyder(a)),theta);

% Using phi and its derivatives, below you should define functions sigma,
% sigmaprime, sigmapprime.

sigma = @(theta) sigma_fun(theta,q_plus,q_tilde,phi);
sigmaprime = @(theta) sigmaprime_fun(theta,q_tilde,phiprime);
sigmapprime = @(theta) sigmapprime_fun(theta,phipprime);

% This is the data structure to be passed to various functions
% You will need to add extra information to it.
data.Kp=Kp;
data.Kd=Kd;
data.D=Dfun;
data.C=Cfun;
data.G=Gfun;
data.B=B;
data.phi=phi;
data.phiprime=phiprime;
data.phipprime=phipprime;
data.sigma=sigma;
data.sigmaprime=sigmaprime;
data.sigmapprime=sigmapprime;

% added info
data.q_tilde = q_tilde;
data.q_plus = q_plus;
data.q_minus = q_minus;

%% HERE WRITE CODE TO TEST WHETHER YOUR VHC WORKS
fprintf('\n Validating VHC...\n')

theta = 0:0.001:1;

% Plot q and check if it passes q plus q minus and q bar
q_test = sigma(theta);
figure
plot(q_test(1,:), q_test(2,:));
xlabel 'q1'
ylabel 'q2'
hold on
title 'VHC q = sigma(theta)'
plot(q_plus(1), q_plus(2),"o");
plot(q_minus(1), q_minus(2),"o");
plot(pi/2,pi,"o");
text(q_plus(1), q_plus(2),"   q plus");
text(q_minus(1), q_minus(2),"   q minus");
text(pi/2,pi,"   q bar");

% exhaustive check that curve is in W
in_W = true;
for i = 2:1000
    q_point = q_test(:,i);

    if q_point(1) > 0 && q_point(1) < pi/2 && ...
       q_point(2) > -2*q_point(1) + 2*pi && q_point(2) < 3*pi
    elseif q_point(1) > pi/2 && q_point(1) < pi && ...
       q_point(2) > -pi && q_point(2) < -2*q_point(1) + 2*pi
    elseif q_point(1) == pi/2 || q_point(1) == 0 || q_point(1) == pi
    else
        in_W = false
    end
end

% verify that the curve is a regular VHC
B_perp = null(B');
regularVHC = [];
for i = 1:1001
    regularVHC = [regularVHC; B_perp'*Dfun(q_test(:,i))*sigmaprime(theta(i))];
end
figure
plot(theta,regularVHC);
ylim([0 0.55])
title 'Regular VHC check'

% verify whether the constraints induces a stable limit cycle
% Building psi
Psi1 = [];
Psi2 = [];
for i = 1:1001
    denominator = B_perp'*Dfun(q_test(:,i))*sigmaprime(theta(i));
    Psi1= [Psi1; -B_perp.'*Gfun(q_test(:,i))/denominator];
    Psi2= [Psi2; -B_perp.'*(Dfun(q_test(:,i))*sigmapprime(theta(i))+Cfun([q_test(:,i);sigmaprime(theta(i))])*sigmaprime(theta(i)))/denominator];
end

% calculate mass, velocity, and delta
M_theta_num = exp(-2*cumtrapz(theta,Psi2));
V_theta_num = -cumtrapz(theta,Psi1.*M_theta_num);

M_theta = @(theta) M_theta_num((theta*1000)+1);
V_theta = @(theta) V_theta_num((theta*1000)+1);

delta_theta = dot(sigmaprime(0), I_delta_num*sigmaprime(1))/...
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

%% NOW YOU CAN SIMULATE THE ROBOT. PLACE YOUR CONTROLLER INSIDE THE FUNCTION acrobot AT THE END OF THIS SCRIPT
ops= odeset('reltol',1e-7,'abstol',1e-7,'Events',@ground_impact);
dt=1/60; % 60 fps; time increment in simulations and animations

fprintf('\n Simulating...\n')
%% DEFINE THE INITIAL CONDITION [q0;qdot0];
% pick initial conditions

simtype = 1;

if simtype == 1 % On the limit cycle
    q0 = sigma(0);
    qdot0 = sigmaprime(0)*delta_theta*sqrt(-2*V_theta(1)/(M_theta(1)-delta_theta^2));

elseif simtype == 2 % On the constraint manifold
    q0 = sigma(0.7);
    qdot0 = sigmaprime(0.7)*delta_theta*sqrt(-2*V_theta(1)/(M_theta(1)-delta_theta^2));

else % outside of constraint manifold
    q0_eps = 0.01;
    qdot0_eps = 0.02;
    q0 = sigma(0.7) + q0_eps;
    qdot0 = sigmaprime(0.7)*delta_theta*sqrt(-2*V_theta(1)/(M_theta(1)-delta_theta^2)) + qdot0_eps;
end

T=[];
X=[];
Te=[];
Ie=[];
Xe=[];
post_impact_state=[q0;qdot0];
% Simulate number_steps steps
for step=1:number_steps
    fprintf('\n...step %d...\n',step);
    [t,x,te,xe,ie]=ode45(@(t,x) acrobot(t,x,data),0:dt:10,post_impact_state,ops);
    % Application of the impact map
    impact_state=xe(end,:)';
    post_impact_state=Deltafun(impact_state);
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
stance_leg=line([ref l*cos(q0(1))],[0 l*sin(q0(1))],'color','red','linewidth',2);
swing_leg=line([ref+l*cos(q0(1)) ref+l*cos(q0(1))+l*cos(q0(1)+q0(2))],...
    [l*sin(q0(1)) l*sin(q0(1))+l*sin(q0(1)+q0(2))],'linewidth',2);
fprintf('\n Animation is ready...\n')


% v = VideoWriter("P1_outside_of_constraints");
% open(v)

animation_slowdown_factor=2; % >1 means slow down
for step=1:length(Ie)
    t=T{step};
    x=X{step};
    xe=Xe{step};
    xe=xe(end,:);
    for k=2:length(t) 
        t0=clock;
        drawnow;
        q=x(k,1:2)';
        xdata1=[ref ref+l*cos(q(1))];
        xdata2=[ref+l*cos(q(1)) ref+l*cos(q(1))+l*cos(q(1)+q(2))];
        ydata1=[0 l*sin(q(1))];
        ydata2=[l*sin(q(1)) l*sin(q(1))+l*sin(q(1)+q(2))];
        set(stance_leg,'xdata',xdata1,'ydata',ydata1);
        set(swing_leg,'xdata',xdata2,'ydata',ydata2);
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
    ref=ref+l*(cos(xe(1))+cos(xe(1)+xe(2)));
end
% close(v);

%% FUNCTIONS
function xdot=acrobot(t,x,data)
q1=x(1);
q2=x(2);
q1dot=x(3);
q2dot=x(4);
q=x(1:2);
qdot=x(3:4);
Kp=data.Kp;
Kd=data.Kd;
D=data.D;
C=data.C;
G=data.G;
B=data.B;
phi=data.phi;
phiprime=data.phiprime;
phipprime=data.phipprime;

q_tilde = data.q_tilde;
q_plus = data.q_plus;
q_minus = data.q_minus;
% DEFINE YOUR CONTROLLER HERE
theta = (q_plus(1) - q1)/q_tilde(1);
theta_dot = -q1dot/q_tilde(1);
e = q2 - phi(theta);
e_dot = q2dot - phiprime(theta)*theta_dot;
tau= inv([phiprime(theta)/q_tilde(1) 1]*inv(D(q))*B)*...
    ([phiprime(theta)/q_tilde(1) 1]*inv(D(q))*(C([q;qdot])*qdot + G(q))+...
    phipprime(theta)*theta_dot^2 - Kp*sin(e) - Kd*e_dot);

qddot=inv(D(q))*(-C(x)*qdot-G(q)+B*tau);
xdot=[qdot;qddot];
end

function [value,isterminal,direction]=ground_impact(t,x)
q1=x(1);
q2=x(2);
% impact occurs when q2 = -2*q1+2*pi
value=q2+2*q1-2*pi;

% We exclude the scuffing point from the impact conditions
if abs(q1-pi/2)<0.01
    isterminal=0;
else
    isterminal=1;
end

% We distinguish between impact on S^+ or S^- by changing the way in which
% ode45 monitors the zero crossing
if abs(q1)<pi/2
    direction=-1;
else
    direction=1;
end
end

% Define sigma functions
function sigma = sigma_fun(theta,q_plus,q_tilde,phi)
    sigma1 = q_plus(1) - theta*q_tilde(1);
    sigma2 = phi(theta);
    sigma = [sigma1; sigma2];
end

function sigmaprime = sigmaprime_fun(theta,q_tilde,phiprime)
    sigma1 = -q_tilde(1);
    sigma2 = phiprime(theta);
    sigmaprime = [sigma1; sigma2];
end

function sigmapprime = sigmapprime_fun(theta,phipprime)
    sigma1 = 0;
    sigma2 = phipprime(theta);
    sigmapprime = [sigma1; sigma2];
end

