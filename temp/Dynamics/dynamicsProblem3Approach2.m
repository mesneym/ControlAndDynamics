%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      PROBLEM 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This approach uses christoffel symbols and accounts for rotational 
% inertia


syms I1_xx I1_xy I1_xz I1_yx I1_yy I1_yz I1_zx I1_zy I1_zz
syms I2_xx I2_xy I2_xz I2_yx I2_yy I2_yz I2_zx I2_zy I2_zz
syms I3_xx I3_xy I3_xz I3_yx I3_yy I3_yz I3_zx I3_zy I3_zz

syms m1 m2 m3 g F1 F2 F3
syms q1(t) q2(t) q3(t) q1_d(t) q1_dd(t) q2_d(t) q2_dd(t) q3_d(t) q3_dd(t)
syms a1c a2c a3c a1 a2 a3 


%%%%%%%%%%%%%%%%%%%%%%
% Symbol definitions
%%%%%%%%%%%%%%%%%%%%%%
% q(t) - Represents generalized coordinates
% F - Represents generalized force or torque
% m - Represents mass of links 
% q_d(t) - Represents derivative of q wrt to t
% q_dd(t)- Represents second derivative of q wrt to t
% r - Position of COM wrt to inertial(world) frame
% a's - robot link parameters
% I - Rotational Inertia from links with respect to inertial frame


I1 = [I1_xx I1_xy I1_xz;
      I1_yx I1_yy I1_yz;
      I1_zx I1_zy I1_zz];

I2 = [I2_xx I2_xy I2_xz;
      I2_yx I2_yy I2_yz;
      I2_zx I2_zy I2_zz];

I3 = [I3_xx I3_xy I3_xz;
      I3_yx I3_yy I3_yz;
      I3_zx I3_zy I3_zz];

Jvc1 =  [-a1c*sin(q1) 0 0;
          a1c*cos(q1) 0 0;
          0           0 0];
 
Jvc2 = [-sin(q1)*(a1+q2+a2c) cos(q1) 0;
         cos(q1)*(a1+q2+a2c) sin(q1) 0;
         0                   0      0];
    
Jvc3 = [-sin(q1)*(a1+q2+a2) cos(q1) -a3c*sin(q3);
         cos(q1)*(a1+q2+a2) sin(q1)  a3c*cos(q3);
         0                   0       0];
    
Jw1 =  [0 0 0;
        0 0 0;
        1 0 0];

Jw2 =  [0 0 0;
        0 0 0;
        0 0 0];
 
Jw3 = [ 0 0 0;
        0 0 0;
        1 0 1];

D =  m1*transpose(Jvc1)*Jvc1 + m2*transpose(Jvc2)*Jvc2 ...
   + m3*transpose(Jvc3)*Jvc3 + transpose(Jw1)*I1*Jw1   ...
   + transpose(Jw2)*I2*Jw2 + transpose(Jw3)*I3*Jw3 ;
D = simplify(D); % inertia matrix

P =     m1*g*a1c*sin(q1) + m2*g*(a1+q2+a2c)*sin(q1) ...
      + m3*g*((a1+q2+a2)*sin(q1) + a3c*sin(q3));

q = [q1; q2; q3];


Dmat = D(t); %This is necessary to index values of D
qmat = q(t); %This is necessary to index values of q

syms C [3 3 3]; %Christofel Symbols is a 3x3x3 matrix

%calculation christofell symbols
for i = 1:3
    for j = 1:3
        for k= 1:3   
              C(i,j,k) = functionalDerivative(Dmat(k,j),qmat(j)) ...
                       + functionalDerivative(Dmat(k,i),qmat(i)) ...     
                       - functionalDerivative(Dmat(i,j),qmat(k));
              C(i,j,k) = sym(1)/2 * C(i,j,k);
        end
    end
end

C = simplify(C);

syms dP [3 1]; %derivate of potential energy wrt q's

for i = 1:3
    dP(i)= functionalDerivative(P,qmat(i));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Printing values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D
C
dP





