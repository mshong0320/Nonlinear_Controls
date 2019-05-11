function f = twolinkarm_feedforward()
clc; clear all; close all;

S.m1 = 1;
S.m2 = 1;
S.l1 = .5;
S.l2 = .5;
S.lc1 = .25;
S.lc2 = .25;
S.I1 = S.m1*S.l1/12;
S.I2 = S.m2*S.l2/12;
S.g = 9.8;

% initial state with high-velocity
x0 = [0; 0; 1; 1; 0; 0];
S.x0 = x0;

% desired state
xd = [5; 7; 0; 0; 0; 0];
S.xd = xd;

% desired accelration -- only applicable for trajectory tracking
S.ad = [0; 0];

T = 20; % simulation time
S.T = T;
S.kp = 20; S.kd = 15;

% perturb initial condition
x0rand = x0+5*rand(size(x0));

[ts, xs] = ode45(@arm_ode, [0 T], x0rand, [], S);

%flat outputs initial and final
y0 = uni_h(x0);
yf = uni_h(xd);

% flat output 1st order
dy0 = [x0(3);x0(4)]; % desired starting velocity
dyf = [xd(3);xd(4)]; % desired end velocity

AA = getA(0, T, x0(1), x0(2), xd(1), xd(2), x0(3), x0(4), xd(3), xd(4));
% AAA = poly3_coeff(y0, dy0, yf, dyf, T);

X = AA*poly3([0:.1:T]);
X_d = AA*dpoly3([0:.1:T]);

%desired acceleration q1_ddot, q2_ddot
X_dd = AA*ddpoly3([0:.1:T]);
Time = [0:0.1:T];

for i = 1:length(X_dd)
    xinput = [X(:,i);
              X_d(:,i)];
    acc = X_dd(:,i);
    
    [M, C, N] = arm_dyn(xinput, S);
    
    b = C*X_d(:,i) + N;
    
    ud(:,i)=M*(acc+inv(M)*b);
    
end


for i = 1:length(ts)
   x = xs(i,1:2)';
   v = xs(i,3:4)';
   a = xs(i,5:6)';
   
   xinput = [x;
             v];
   
   [M, C, N] = arm_dyn(xinput, S);
   b = C*v + N;
   u_exec(:,i) = M*(a+inv(M)*b);
end

figure(1)
plot(X(1,:),X(2,:), '-r')

hold on
plot(xs(:,1), xs(:,2), '-b')
plot(S.xd(1),S.xd(2),'*g')
xlabel('shoulder q1')
ylabel('elbow q2')
legend('qd(t)', 'executed', 'goal')
hold off

figure(2)
plot(Time, ud(1,:), '-r')
hold on
plot(Time, ud(2,:), '-b')
plot(ts, u_exec(1,:), '-m')
plot(ts, u_exec(2,:), '-k')
xlabel('time')
ylabel('Joint angles in radians')
legend('ud1(t)', 'ud2(t)', 'u1(t)', 'u2(t)')

% for i = 1:length(ts)
%     t = ts(i);
%     actual_xdd = AA*ddpoly3(t);
%     u_d = S.M*(actual_xdd)
% %     u = 
% end


function dx = arm_ode(t, x, S)
% the ODE for the arm
Y = [0 S.x0(3) S.xd(1) 0;
     0 S.x0(4) S.xd(2) 0];

LambdaInv = [1 0 -3/S.T^2 2/S.T^3;
             0 1 -2/S.T 1/S.T^2;
             0 0 3/S.T^2 -2/S.T^3;
             0 0 -1/S.T 1/S.T^2];

A_mat = Y*LambdaInv;
lambda = [ones(size(t));
          t;
          t.^2;
          t.^3];

y = A_mat * lambda;

dpoly3 = [zeros(size(t));
          ones(size(t));
          2*t;
          3*t.^2];
      
ddpoly3 = [zeros(size(t));
           zeros(size(t));
           2;
           6*t];

y_dot = A_mat*dpoly3;
y_ddot = A_mat*ddpoly3;

[M, C, N] = arm_dyn(x, S);

q = x(1:2);
v = x(3:4);

b = C*v + N;
virtualinput = y_ddot-S.kp*(q-y)-S.kd*(v-y_dot);

u = M*virtualinput + b;
acc = -inv(M)*b+inv(M)*u;

dx = [v;
      inv(M)*(u - C*v - N);
      acc];


function [M, C, N] = arm_dyn(x, S)
% compute the dynamical terms of the manipulator
% M - mass matrix
% C - Corilois matrix
% N - gravity/damping forces vector
q = x(1:2);
v = x(3:4);

c1 = cos(q(1));
c2 = cos(q(2));
s2 = sin(q(2));
c12 = cos(q(1) + q(2));

% coriolis matrix
C = -S.m2*S.l1*S.lc2*s2*[v(2), v(1) + v(2);
                        -v(1), 0] + diag([.2;.2]);
% mass elements
m11 = S.m1*S.lc1^2 + S.m2*(S.l1^2 + S.lc2^2 + 2*S.l1*S.lc2*c2) + ...
      S.I1 + S.I2;
m12 = S.m2*(S.lc2^2 + S.l1*S.lc2*c2) + S.I2;
m22 = S.m2*S.lc2^2 + S.I2;

% mass matrix
M = [m11, m12;
     m12, m22];

% gravity, damping, etc...
N = [(S.m1*S.lc1 + S.m2*S.l1)*S.g*c1 + S.m2*S.lc2*S.g*c12;
      S.m2*S.lc2*S.g*c12];
  
  
function A = getA(t0, tf, p0x, p0y, pfx, pfy, v0x, v0y, vfx, vfy)
Y = [p0x v0x pfx vfx;
     p0y v0y pfy vfy];

lambda_i = [1; t0; t0^2; t0^3];
lambda_f = [1; tf; tf^2; tf^3];

lambda_i_dot = [0; 1; 2*t0; 3*t0^2];
lambda_f_dot = [0; 1; 2*tf; 3*tf^2];

L = [lambda_i lambda_i_dot lambda_f lambda_f_dot];
A = Y*inv(L);


function A = poly3_coeff(y0, dy0, yf, dyf, T)
% computes cubic curve connecting (y0,dy0) and (yf, dyf) at time T
Y = [y0, dy0, yf, dyf];
L = [poly3(0), dpoly3(0), poly3(T), dpoly3(T)];
A = Y*inv(L);

function f = poly3(t)
f = [ones(size(t)); t; t.^2; t.^3];

function f = dpoly3(t)
f = [zeros(size(t)); ones(size(t)); 2*t; 3*t.^2];

function f= ddpoly3(t)
f = [zeros(size(t)); zeros(size(t)); 2*ones(size(t)); 6*t];

function y = uni_h(x)
% output function
y = x(1:2);

