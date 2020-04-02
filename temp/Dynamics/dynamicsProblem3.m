%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      PROBLEM 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note This approach doesn't account for Rotational Inertia'

syms I  m1 m2 m3 g F1 F2 F3
syms q1(t) q2(t) q3(t) q1_d(t) q1_dd(t) q2_d(t) q2_dd(t) q3_d(t) q3_dd(t)
syms a1c a2c a3c a1 a2 a3 

%%%%%%%%%%%%%%%%%%%%%%
% Symbol definitions
%%%%%%%%%%%%%%%%%%%%%%
% q(t) - Represents generalized coordinates
% F - Represents generalized force or torque
% q_d(t) - Represents derivative of q wrt to t
% q_dd(t)- Represents second derivative of q wrt to t
% r - Position of COM wrt to inertial(world) frame
% a's - robot link parameters
% I - Rotational Inertia from links

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Lagrangian 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r1 = [a1c*cos(q1)
      a1c*sin(q1)];

r2 = [(a1+q2+a2c)*cos(q1)
      (a1+q2+a2c)*sin(q1)];

r3 = [(a1+ q2 + a2)*cos(q1) + a3c*cos(q3)
      (a1+ q2 + a2)*sin(q1) + a3c*sin(q3)];
  
% Getting velocities
v1 = diff(r1);  
v2 = diff(r2);
v3 = diff(r3);

% Kinetic Energy
T1 = sym(1)/2 * transpose(v1) * v1;
T2 = sym(2)/2 * transpose(v2) * v2;
T3 = sym(3)/3 * transpose(v3) * v3;

T = T1 + T2 + T3

% Potential Energy
V1 = m1*g*a1c*sin(q1);
V2 = m2*g*(a1c + q2* + a2c)*sin(q1);
V3 = m3*g*((a1+q2*+a2)*sin(q1) + a3c*sin(q3));

V = V1 + V2 + V3

% Lagrangian
L = T-V


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Taking derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Derivative of L wrt to q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dL_dq1 = functionalDerivative(L,q1);
dL_dq2 = functionalDerivative(L,q2);
dL_dq3 = functionalDerivative(L,q3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Derivative of L wrt q_dot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Replacing derivatives with symbols
% Matlab can't take functional derivative with
% respect to a derivative 
L=subs(L,diff(q1(t),t),q1_d);
L=subs(L,diff(q2(t),t),q2_d);
L=subs(L,diff(q3(t),t),q3_d);

dL_dq1_dot= functionalDerivative(L,q1_d);
dL_dq2_dot= functionalDerivative(L,q2_d);
dL_dq3_dot= functionalDerivative(L,q3_d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Derivative of dl/dq_dot wrt to t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Taking derivative of dL_dq with respect to t
dt_dL_dq1_dot = diff(dL_dq1_dot,t);
dt_dL_dq2_dot = diff(dL_dq2_dot,t);
dt_dL_dq3_dot = diff(dL_dq3_dot,t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Equation of Motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replacing derivatives with symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Replacing derivatives with symbols to perform substaction
dt_dL_dq1_dot=subs(dt_dL_dq1_dot,diff(q1_d(t),t),q1_dd);
dt_dL_dq1_dot=subs(dt_dL_dq1_dot,diff(q1(t),t),q1_d);
dt_dL_dq1_dot=subs(dt_dL_dq1_dot,diff(q2_d(t),t),q2_dd);
dt_dL_dq1_dot=subs(dt_dL_dq1_dot,diff(q2(t),t),q2_d);
dt_dL_dq1_dot=subs(dt_dL_dq1_dot,diff(q3_d(t),t),q3_dd);
dt_dL_dq1_dot=subs(dt_dL_dq1_dot,diff(q3(t),t),q3_d);

dt_dL_dq2_dot=subs(dt_dL_dq2_dot,diff(q1_d(t),t),q1_dd);
dt_dL_dq2_dot=subs(dt_dL_dq2_dot,diff(q1(t),t),q1_d);
dt_dL_dq2_dot=subs(dt_dL_dq2_dot,diff(q2_d(t),t),q2_dd);
dt_dL_dq2_dot=subs(dt_dL_dq2_dot,diff(q2(t),t),q2_d);
dt_dL_dq2_dot=subs(dt_dL_dq2_dot,diff(q3_d(t),t),q3_dd);
dt_dL_dq2_dot=subs(dt_dL_dq2_dot,diff(q3(t),t),q3_d);

dt_dL_dq3_dot=subs(dt_dL_dq3_dot,diff(q1_d(t),t),q1_dd);
dt_dL_dq3_dot=subs(dt_dL_dq3_dot,diff(q1(t),t),q1_d);
dt_dL_dq3_dot=subs(dt_dL_dq3_dot,diff(q2_d(t),t),q2_dd);
dt_dL_dq3_dot=subs(dt_dL_dq3_dot,diff(q2(t),t),q2_d);
dt_dL_dq3_dot=subs(dt_dL_dq3_dot,diff(q3_d(t),t),q3_dd);
dt_dL_dq3_dot=subs(dt_dL_dq3_dot,diff(q3(t),t),q3_d);


% Replacing derivatives with symbols to perform substaction
dL_dq1=subs(dL_dq1,diff(q1(t),t,t),q1_dd);
dL_dq1=subs(dL_dq1,diff(q1(t),t),q1_d);
dL_dq1=subs(dL_dq1,diff(q2(t),t,t),q2_dd);
dL_dq1=subs(dL_dq1,diff(q2(t),t),q2_d);
dL_dq1=subs(dL_dq1,diff(q3(t),t,t),q3_dd);
dL_dq1=subs(dL_dq1,diff(q3(t),t),q3_d);

dL_dq2=subs(dL_dq2,diff(q1(t),t,t),q1_dd);
dL_dq2=subs(dL_dq2,diff(q1(t),t),q1_d);
dL_dq2=subs(dL_dq2,diff(q2(t),t,t),q2_dd);
dL_dq2=subs(dL_dq2,diff(q2(t),t),q2_d);
dL_dq2=subs(dL_dq2,diff(q3(t),t,t),q3_dd);
dL_dq2=subs(dL_dq2,diff(q3(t),t),q3_d);

dL_dq3=subs(dL_dq3,diff(q1(t),t,t),q1_dd);
dL_dq3=subs(dL_dq3,diff(q1(t),t),q1_d);
dL_dq3=subs(dL_dq3,diff(q2(t),t,t),q2_dd);
dL_dq3=subs(dL_dq3,diff(q2(t),t),q2_d);
dL_dq3=subs(dL_dq3,diff(q3(t),t,t),q3_dd);
dL_dq3=subs(dL_dq3,diff(q3(t),t),q3_d);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %        Finding Equations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eqn1 = simplify(dt_dL_dq1_dot - dL_dq1) == F1
eqn2 = simplify(dt_dL_dq2_dot - dL_dq2) == F2
eqn3 = simplify(dt_dL_dq3_dot - dL_dq3) == F3
