function q = ik_leg(pE, P, branch1, branch3)
%IK_LEG Closed-form inverse kinematics of the 3-DoF leg.
%   q = IK_LEG(pE, P) returns q = [q1 q2 q3] (rad) that place the foot at
%   pE = [px; py; pz] (cm). branch1 = +1/-1 picks the solution for q1 and
%   branch3 = +1/-1 the sign of cos(q3); the defaults (+1, +1) give the
%   configuration q = 0 at the stance point A of the default task.
%
%   1) l1 = c1*px + s1*py. With t = tan(q1/2):
%        q1 = 2*atan2(py + branch1*sqrt(px^2 + py^2 - l1^2), l1 + px)
%   2) K = s1*px - c1*py = l3*c23 - l2*s2,   H = pz - l0 = l2*c2 + l3*s23
%        s3 = (K^2 + H^2 - l2^2 - l3^2) / (2*l2*l3),  c3 = branch3*sqrt(1 - s3^2)
%   3) with a = l3*c3 and b = l3*s3 + l2:  K = a*c2 - b*s2,  H = b*c2 + a*s2
%        c2 = (a*K + b*H)/(a^2 + b^2),   s2 = (a*H - b*K)/(a^2 + b^2)

if nargin < 3, branch1 = 1; end
if nargin < 4, branch3 = 1; end
px = pE(1); py = pE(2); pz = pE(3);

r2 = px^2 + py^2 - P.l1^2;
if r2 < 0
    error('ik_leg:unreachable', 'Point (%.2f, %.2f, %.2f) is closer than l1 to the joint-1 axis.', px, py, pz);
end
q1 = 2*atan2(py + branch1*sqrt(r2), P.l1 + px);

K = sin(q1)*px - cos(q1)*py;
H = pz - P.l0;
s3 = (K^2 + H^2 - P.l2^2 - P.l3^2) / (2*P.l2*P.l3);
if abs(s3) > 1
    error('ik_leg:unreachable', 'Point (%.2f, %.2f, %.2f) is outside the workspace.', px, py, pz);
end
c3 = branch3*sqrt(1 - s3^2);
q3 = atan2(s3, c3);

a = P.l3*c3;  b = P.l3*s3 + P.l2;
q2 = atan2(a*H - b*K, a*K + b*H);   % the common factor 1/(a^2+b^2) > 0 cancels

q = [q1, q2, q3];
end
