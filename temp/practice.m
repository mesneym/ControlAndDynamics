%Matlab Code
syms alpha theta a d

A(theta,d,a,alpha)= [cos(theta), -sin(theta)*cos(alpha) ...
                     sin(theta)*sin(alpha), a*cos(theta);
                     sin(theta),  cos(theta)*cos(alpha), ... 
                     -cos(theta)*sin(alpha), a*sin(theta);
                     0 sin(alpha) cos(alpha) d;
                     0 0 0 1];
       
syms d1 theta1 theta2 theta3 theta4 l2 l4 l5
A_0_1 = simplify(A(pi/2,d1,0,pi/2));
A_1_2 = simplify(A(theta1+pi/2, l2, 0, pi/2));
A_2_3 = simplify(A(theta2+pi/2, 0 , 0 ,pi/2));
A_3_4 = simplify(A(theta3+pi/2, l4, 0 ,pi/2));
A_4_5 = simplify(A(theta4+pi/2, 0, 0, pi/2));
A_5_n = simplify(A(0, 0, l5, 0));

T_O_n=A_0_1*A_1_2*A_2_3*A_3_4*A_4_5*A_5_n;

simplify(T_O_n)



