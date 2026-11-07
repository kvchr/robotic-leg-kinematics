function [pE, joints, T0E, frames] = fk_leg(q, P)
%FK_LEG Forward kinematics of the 3-DoF leg (Denavit-Hartenberg model).
%   [pE, joints, T0E] = FK_LEG(q, P) for joint angles q = [q1 q2 q3] (rad)
%   returns the foot position pE (3x1, cm) in the base frame {0}, the
%   positions of the points used to draw the leg (3x5: O0, joint-1 frame
%   origin, joint-2 axis, knee, foot), the homogeneous transform T0E and
%   the transforms of all D-H frames {0}, {1}, {2}, {3}, {E}.
%
%   D-H table (theta, d, a, alpha):
%     1: q1 + pi/2,  l0,  0,  -pi/2
%     2: q2 - pi/2, -l1,  l2,  0
%     3: q3,          0,  0,  -pi/2
%     E: pi/2,      -l3,  0,   pi/2

q1 = q(1); q2 = q(2); q3 = q(3);
A1 = dh(q1 + pi/2,  P.l0, 0,    -pi/2);
A2 = dh(q2 - pi/2, -P.l1, P.l2,  0);
A3 = dh(q3,         0,    0,    -pi/2);
AE = dh(pi/2,      -P.l3, 0,     pi/2);

T01 = A1;  T02 = T01*A2;  T03 = T02*A3;  T0E = T03*AE;
pE = T0E(1:3, 4);

% Point on the joint-2 axis at the end of the offset l1 (drawn as a link).
pJ2 = T01(1:3, 4) - P.l1*T01(1:3, 3);
joints = [zeros(3,1), T01(1:3,4), pJ2, T02(1:3,4), pE];
frames = {eye(4), T01, T02, T03, T0E};   % frames {0}, {1}, {2}, {3}, {E}
end

function A = dh(theta, d, a, alpha)
ct = cos(theta); st = sin(theta); ca = cos(alpha); sa = sin(alpha);
A = [ct, -st*ca,  st*sa, a*ct;
     st,  ct*ca, -ct*sa, a*st;
     0,   sa,     ca,    d;
     0,   0,      0,     1];
end
