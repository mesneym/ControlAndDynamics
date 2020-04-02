

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      PROBLEM 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


syms l1 l2 q1(t) q2(t) r M g F1 F2  q1_d(t) q1_dd(t) q2_d(t) q2_dd(t)

%%%%%%%%%%%%%%%%%%%%%%
% Symbol definitions
%%%%%%%%%%%%%%%%%%%%%%
% q(t) - Represents theta(t)
% F - Represents generalized force or torque
% q_d(t) - Represents derivative of q wrt to t
% q_dd(t)- Represents second derivative of q wrt to t
% r - Position of COM wrt to inertial(world) frame
% l - link length of robot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Lagrangian 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = [(l2*cos(q2)+l1)*cos(q1)
     (l2*cos(q2)+l1)*sin(q1)
     l2*sin(q2)]

v = diff(r);  %derivate of r
T = simplify(sym(1)/2* M * transpose(v)*v)
V = M*g*l2*sin(q2)
L = T-V



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Taking derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Derivative of L wrt to q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dL_dq1 = functionalDerivative(L,q1);
dL_dq2 = functionalDerivative(L,q2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Derivative of L wrt q_dot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Replacing derivatives with symbols
% Matlab can't take functional derivative with
% respect to a derivative 
L=subs(L,diff(q2(t),t),q2_d);
L=subs(L,diff(q1(t),t),q1_d);

dL_dq1_dot= functionalDerivative(L,q1_d);
dL_dq2_dot= functionalDerivative(L,q2_d);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Derivative of dl/dq_dot wrt to t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Taking derivative of dL_dq with respect to t
dt_dL_dq1_dot = diff(dL_dq1_dot,t);
dt_dL_dq2_dot = diff(dL_dq2_dot,t);



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


dt_dL_dq2_dot=subs(dt_dL_dq2_dot,diff(q1_d(t),t),q1_dd);
dt_dL_dq2_dot=subs(dt_dL_dq2_dot,diff(q1(t),t),q1_d);
dt_dL_dq2_dot=subs(dt_dL_dq2_dot,diff(q2_d(t),t),q2_dd);
dt_dL_dq2_dot=subs(dt_dL_dq2_dot,diff(q2(t),t),q2_d);



% Replacing derivatives with symbols to perform substaction
dL_dq1=subs(dL_dq1,diff(q1(t),t,t),q1_dd);
dL_dq1=subs(dL_dq1,diff(q1(t),t),q1_d);
dL_dq1=subs(dL_dq1,diff(q2(t),t,t),q2_dd);
dL_dq1=subs(dL_dq1,diff(q2(t),t),q2_d);


dL_dq2=subs(dL_dq2,diff(q1(t),t,t),q1_dd);
dL_dq2=subs(dL_dq2,diff(q1(t),t),q1_d);
dL_dq2=subs(dL_dq2,diff(q2(t),t,t),q2_dd);
dL_dq2=subs(dL_dq2,diff(q2(t),t),q2_d);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Finding Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eqn1 = simplify(dt_dL_dq1_dot - dL_dq1);
enq1=simplify(subs(eqn1,l1,2*l2))== F1

eqn2 =simplify(dt_dL_dq2_dot - dL_dq2);
eqn2 = simplify(subs(eqn2,l1,2*l2))== F2



