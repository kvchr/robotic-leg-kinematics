function [JL, detJ] = jacobian_leg(q, P)
%JACOBIAN_LEG Linear-velocity Jacobian of the foot, v = JL * qdot.
%   [JL, detJ] = JACOBIAN_LEG(q, P) returns the 3x3 Jacobian and its
%   determinant  det(JL) = l2*l3*cos(q3)*(l3*cos(q2+q3) - l2*sin(q2)).

s1 = sin(q(1)); c1 = cos(q(1));
s2 = sin(q(2)); c2 = cos(q(2));
s23 = sin(q(2)+q(3)); c23 = cos(q(2)+q(3));
l1 = P.l1; l2 = P.l2; l3 = P.l3;

% Foot offset from the point (c1*l1, s1*l1, l0) of the joint-2 axis:
% (s1*K, -c1*K, H), i.e. K along (s1, -c1, 0), perpendicular to z0, and H along z0.
K = l3*c23 - l2*s2;
H = l2*c2 + l3*s23;   % = pz - l0

JL = [-s1*l1 + c1*K,  -s1*H,  -s1*l3*s23;
       c1*l1 + s1*K,   c1*H,   c1*l3*s23;
       0,              K,      l3*c23];
detJ = l2*l3*cos(q(3))*K;
end
