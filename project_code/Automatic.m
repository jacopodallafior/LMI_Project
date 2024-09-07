%% INITIALIZATION
clear
close all
clc

%% Parameters
mp = 0.1;   % mass of the pendulum [kg]
Lp = 0.3;   % length of the pendulum [m]
lp = 0.15;  % position of the center of mass of the pendulum [m]
Jp = 0.003; % moment of inertia of the pendulum [kg*m^2]
Lr = 0.15;  % length of the rotary arm [m]
Jr = 0.01;  % moment of inertia of the rotary arm [kg*m^2]
Bp = 0.4;   % viscous friction coefficient of the pendulum [Nms/rad]
Br = 0.2;   % viscous friction coefficient of the rotary arm [Nms/rad]
g = 9.81;   % gravity constant [m/s^2]
%%
% -- Simulation settings
simT = 0:0.01:6; %integrate over this time Time Domain
q0   = [0.05; 0; 0.06; 0]; %Initial Conditions

%% Equilibrium point

theta_eq = 0;
alpha_eq = 0;
theta_dot_eq = 0;
alpha_dot_eq = 0;

% Define symbolic variables
syms theta alpha theta_dot alpha_dot u

% Define the state vector and input
x = [theta; alpha; theta_dot; alpha_dot];


% M(q)
M = [Jr + mp*(Lr^2 + lp^2*(1-cos(alpha))^2), mp*lp*Lr*cos(alpha);
    mp*lp*Lr*cos(alpha), Jp + mp*lp^2];

% C(q, q_dot)
C = [2*mp*lp^2*alpha_dot*sin(alpha)*cos(alpha), -mp*lp*Lr*alpha_dot*sin(alpha);
    -mp*lp^2*theta_dot*sin(alpha)*cos(alpha), 0];
disp(' M = ');
disp(M);
disp(' C = ');
disp(C);

fv = [Br*theta_dot;
    Bp*alpha_dot];

% G(q)
G = [0;
    -mp*lp*g*sin(alpha)];

% tau
tau = [u; 0];

% Define the equations of motion
q = [theta; alpha];
q_dot = [theta_dot; alpha_dot];
q_ddot = M \ (tau - C*q_dot - fv - G);

% State-space 
f = [q_dot;
    q_ddot];

% Linearization at equilibria

A = double(subs(jacobian(f, x), {theta, alpha, theta_dot, alpha_dot, u}, {theta_eq, alpha_eq, theta_dot_eq, alpha_dot_eq, 0}));
B = double(subs(jacobian(f, u), {theta, alpha, theta_dot, alpha_dot, u}, {theta_eq, alpha_eq, theta_dot_eq, alpha_dot_eq, 0}));

% Display the results
disp('Linearized state-space:');
disp('A = ');
disp(A);
disp('B = ');
disp(B);

% -- Plant dimensions
[n,~] = size(A); %[n x n]
[~,p] = size(B); %[n x p]

%% Stability: Eigenvalues Test

% Eigenvalue to check stability

ev = real(eig(A)); %real part of the eigenvalues of A

disp('--------------------------------------------------------------------------')
disp('Internal Stability Test (open-loop)')
disp('*Eigenvalues test:')
if all(ev<0)
    disp(['    ' 'all eigenvalues have negative real part --> system is stable'])
else
    disp(['    ' 'at least one eigenvalue has positive real part --> system is unstable'])
end
disp('*Open-loop eigenvalues real parts: ')
disp(['    ', mat2str(ev')]);

%% Controllability Test

% Check the controllability using C matrix

ctr = rank(ctrb(A,B)); %controllability matrix
Control=ctrb(A,B);
disp('--------------------------------------------------------------------------')
disp('Controllability Test');
disp('*Controllability matrix test:');
if ctr==n
    disp(['    ' 'ctrb matrix is full rank --> system is controllable'])
else
    disp(['    ' 'ctrb matrix is rank deficient --> system is NOT controllable'])
    disp('*plant dimension:');
    disp(['    ', num2str(n)]);
    disp('*ctrb matrix rank:');
    disp(['    ', num2str(ctr), newline]);
end


%% Feedback design: desired Convergence Rate

% -- Clear the internal memory of YALMIP
yalmip('clear')

% -- Define settings to give to the solver
opts = sdpsettings('solver', 'mosek', 'verbose', 0);

fake_zero = 1e-3;

% -- Optimization variables
W = sdpvar(n,n);
X = sdpvar(p,n);

% -- Constraints
alpha_b = 2;
constr = [ W >= fake_zero*eye(n)                ;...
    (A*W+B*X)+(A*W+B*X)' <= -2*alpha_b*W ];

% -- Solve the problem
sol = solvesdp(constr,[],opts);

% -- extract the results
W_b = value(W);
X_b = value(X);
K1 = X_b*inv(W_b);

% -- Check constraints
[primalfeas,dualfeas] = check(constr);
disp('--------------------------------------------------------------------------')
disp(['LMI Feedback design (alpha = ' num2str(alpha_b) ')'])
disp('*Constraints check:')
if all([primalfeas;dualfeas] > -fake_zero)
    disp(['    ', 'Constraints are satisfied --> system is stable'])
    disp('*Gain matrix K:');
    disp(['    ', num2str(K1)]);
else
    disp(['    ', 'Constraints are NOT satisfied --> system is unstable'])
end
disp('*Closed-loop eigenvalues real parts:')
disp(['    ' num2str(real(eig(A+B*K1)'))]);

% -- Simulation
[simT_b,simQ_b]  = ode45(@(t,q)f_NLDyna(q, Jr, mp, Lr, lp, Jp, Br, Bp, g, K1*q),simT,q0);
% M value of first control
P_b = inv(W_b);
M_b = sqrt(max(eig(P_b))/min(eig(P_b)));

%% Feedback design: desired Convergence Rate and minimum control effort
% -- Clear the internal memory of YALMIP
yalmip('clear')

% -- Optimization variables
W = sdpvar(n,n);
X = sdpvar(p,n);
k = sdpvar(1,1);

% -- Constraints
alpha_c = 2;
constr = [ W >= eye(n)               ;...
    (A*W+B*X)+(A*W+B*X)' <= -2*alpha_c*W ;...
    [ k*eye(n)      X'   ;...
    X      k*eye(p) ]   >= fake_zero*eye(n+p)];

% -- Solve the problem
sol = solvesdp(constr,k,opts);

% -- Extract the results
W_c = value(W);
X_c = value(X);
k_c = value(k);
K2 = X_c*inv(W_c);

% -- Check constraints
[primalfeas,dualfeas] = check(constr);
disp('--------------------------------------------------------------------------')
disp(['LMI Feedback design (alpha = ' num2str(alpha_c) ', minimize |K|)' ])
disp('*Constraints check:')
if all([primalfeas;dualfeas] > -fake_zero)
    disp(['    ' 'Constraints are satisfied --> system is stable'])
    disp('*Gain matrix K:');
    disp(['    ' num2str(K2)]);
    disp(['*Minimum bound on gain matrix norm: |K| < ', num2str(k_c)]);
else
    disp(['    ', 'Constraints are NOT satisfied --> system is unstable'])
end
disp('Closed-loop eigenvalues real parts:')
disp(['    ', num2str(real(eig(A+B*K2)'))]);

% -- Simulation
[simT_c,simQ_c]  = ode45(@(t,q)f_NLDyna(q, Jr, mp, Lr, lp, Jp, Br, Bp, g, K2*q),simT,q0);
% M value of second control
M_c = sqrt(k_c);


%% Estimation of the overshoot M1

yalmip('clear');

alpha_bar = 2;

% -- variables
P = sdpvar(n,n);
k = sdpvar(1,1);

% -- LMIs
constr = [ eye(n) <= P <= k*eye(n);
    (A+B*K1)'*P + P*(A+B*K1) <= -2*alpha_bar*P ];

sol = solvesdp(constr, k, opts);

[primal_res, dual_res] = check(constr);

feas = all([primal_res; dual_res] > -fake_zero);


disp(newline);
disp('-----------------------')
disp('Estimating overshoot M');

k_val = value(k);
M1 = sqrt(k_val);

if(feas)
    disp( ['The problem is feasible, value of M: ', mat2str(M1)] );
    disp('-----------------------')
else
    disp('The problem is infeasible');
    disp('-----------------------')
end


%% Estimation of the overshoot M2

yalmip('clear');

% -- variables
P = sdpvar(n,n);
k = sdpvar(1,1);

% -- LMIs
constr = [ eye(n) <= P <= k*eye(n);
    (A+B*K2)'*P + P*(A+B*K2) <= -2*alpha_bar*P ];

sol = solvesdp(constr, k, opts);

[primal_res, dual_res] = check(constr);

feas = all([primal_res; dual_res] > -fake_zero);


disp(newline);
disp('-----------------------')
disp('Estimating overshoot M');

k_val = value(k);
M2 = sqrt(k_val);

if(feas)
    disp( ['The problem is feasible, value of M: ', mat2str(M2)] );
else
    disp('The problem is infeasible');
    disp('-----------------------')
end

%% Compare solutions
% -- State norm evolution over time

normX_b = vecnorm(simQ_b');
normX_c = vecnorm(simQ_c');
% Exponential bounds
x_bound_b = M_b*norm(q0)*exp(-alpha_bar*simT);
x_bound_c = M_c*norm(q0)*exp(-alpha_bar*simT);

% -- Input evolution over time

u_b = (K1*simQ_b')';
u_c = (K2*simQ_c')';


%% Plot 1
figure(1)


plot(simT_b,normX_b, 'Color',[0 .4470 .7410], 'LineWidth',1.5)
hold on
plot(simT_c,normX_c, 'Color',[.9290 .6940 .1250], 'LineWidth',1.5)
grid on
grid minor
legend('$\alpha$ only','minimize $|K|$', 'interpreter', 'latex')
title('State Norm')
xlabel('$t[s]$', 'interpreter', 'latex')

figure(2), clf
hold on
grid minor
plot(simT_b,normX_b, 'Color',[0 .4470 .7410], 'LineWidth',1.5)
plot(simT_b, x_bound_b, '--r');
legend('$|x(t)|$', '$M|x(0)|e^{-\bar \alpha t}$', 'interpreter', 'latex');

xlabel('$t[s]$', 'interpreter', 'latex');

figure(3)
hold on
grid minor
plot(simT_c,normX_c, 'Color',[0 .4470 .7410], 'LineWidth',1.5)
plot(simT_c, x_bound_c, '--r');
legend('$|x(t)|$', '$M|x(0)|e^{-\bar \alpha t}$', 'interpreter', 'latex');

xlabel('$t[s]$', 'interpreter', 'latex');


%% Plot 2
figure(4)

plot(simT_b,simQ_b(:,1), 'Color',[0 .4470 .7410], 'LineWidth',1.5)
hold on
plot(simT_c,simQ_c(:,1), 'Color',[.9290 .6940 .1250], 'LineWidth',1.5)
grid on
grid minor
legend('$\alpha$ only','minimize $|K|$', 'interpreter', 'latex')
title('theta Position')
xlabel('$t[s]$', 'interpreter', 'latex')

figure(5)

plot(simT_b,simQ_b(:,2), 'Color',[0 .4470 .7410], 'LineWidth',1.5)
hold on
plot(simT_c,simQ_c(:,2), 'Color',[.9290 .6940 .1250], 'LineWidth',1.5)
grid on
grid minor
legend('$\alpha$ only','minimize $|K|$', 'interpreter', 'latex')
title('alpha position')
xlabel('$t[s]$', 'interpreter', 'latex')

%% Plot 3

figure(6)

plot(simT_b,u_b, 'Color',[0 .4470 .7410], 'LineWidth',1.5)
hold on
plot(simT_c,u_c, 'Color',[.9290 .6940 .1250], 'LineWidth',1.5)
grid on
grid minor
legend('$\alpha$ only','minimize $|K|$', 'interpreter', 'latex')
title('Control Input')
xlabel('$t[s]$', 'interpreter', 'latex')
