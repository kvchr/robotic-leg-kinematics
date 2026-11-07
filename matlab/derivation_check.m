%DERIVATION_CHECK Symbolic check of the kinematic model (Part A).
%   Builds the D-H transforms, prints the forward kinematics, the Jacobian
%   and its determinant, and compares them with the hand-written formulas
%   used in fk_leg.m, jacobian_leg.m and ik_leg.m. It then tests the
%   closed-form inverse kinematics on random configurations.
%   Needs the Symbolic Math Toolbox. Output: ../results/derivation_check.txt

here = fileparts(mfilename('fullpath'));
outDir = fullfile(here, '..', 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
logFile = fullfile(outDir, 'derivation_check.txt');
if exist(logFile, 'file'), delete(logFile); end
diary(logFile);


syms q1 q2 q3 l0 l1 l2 l3 real
dh = @(th,d,a,al) [cos(th) -sin(th)*cos(al)  sin(th)*sin(al) a*cos(th);
                   sin(th)  cos(th)*cos(al) -cos(th)*sin(al) a*sin(th);
                   0        sin(al)          cos(al)         d;
                   0 0 0 1];
A1 = dh(sym(pi)/2 + q1,  l0, 0,  -sym(pi)/2);
A2 = dh(q2 - sym(pi)/2, -l1, l2,  0);
A3 = dh(q3,              0,  0,  -sym(pi)/2);
AE = dh(sym(pi)/2,      -l3, 0,   sym(pi)/2);
T0E = simplify(A1*A2*A3*AE);

fprintf('Rotation R_E^0:\n'); disp(T0E(1:3,1:3));
p = T0E(1:3,4);
fprintf('Foot position p_E:\n'); disp(p);

K = l3*cos(q2+q3) - l2*sin(q2);
H = l2*cos(q2) + l3*sin(q2+q3);
p_formula = [cos(q1)*l1 + sin(q1)*K; sin(q1)*l1 - cos(q1)*K; l0 + H];
fprintf('p_E minus [c1*l1 + s1*K; s1*l1 - c1*K; l0 + H]:\n'); disp(simplify(p - p_formula));

J = jacobian(p, [q1 q2 q3]);
J_formula = [-sin(q1)*l1 + cos(q1)*K, -sin(q1)*H, -sin(q1)*l3*sin(q2+q3);
              cos(q1)*l1 + sin(q1)*K,  cos(q1)*H,  cos(q1)*l3*sin(q2+q3);
              0,                       K,          l3*cos(q2+q3)];
fprintf('J_L minus the formula in jacobian_leg.m:\n'); disp(simplify(J - J_formula));

% Angular part: the joint axes z0, z1, z2 in the base frame.
T01 = A1; T02 = simplify(A1*A2);
fprintf('Joint axes [z0 z1 z2]:\n'); disp([sym([0;0;1]) T01(1:3,3) T02(1:3,3)]);

d = simplify(det(J));
fprintf('det(J_L) factors:\n'); disp(factor(d));
fprintf('det(J_L) - l2*l3*cos(q3)*K = %s\n', char(simplify(d - l2*l3*cos(q3)*K)));

% Closed-form inverse of J_L (used in the report).
syms vx vy vz real
c1 = cos(q1); s1 = sin(q1); c3 = cos(q3); s23 = sin(q2+q3); c23 = cos(q2+q3);
w = -s1*vx + c1*vy;
qd1 = (c1*vx + s1*vy)/K;
qd2 = (c23*(w - l1*qd1) - s23*vz)/(l2*c3);
qd3 = (H*vz - K*(w - l1*qd1))/(l2*l3*c3);
fprintf('J_L*[qd1; qd2; qd3] - v (closed-form inverse):\n');
disp(simplify(J*[qd1; qd2; qd3] - [vx; vy; vz]));

% Numeric test of the closed-form inverse kinematics: every one of the
% four branches must place the foot back at the target.
P = leg_params();
rng(0);
nTest = 5000; maxErr = 0;
for k = 1:nTest
    q = (rand(1,3) - 0.5)*2*pi;
    pE = fk_leg(q, P);
    for b1 = [1 -1]
        for b3 = [1 -1]
            qs = ik_leg(pE, P, b1, b3);
            maxErr = max(maxErr, norm(fk_leg(qs, P) - pE));
        end
    end
end
fprintf(['Inverse kinematics: max position error over %d random configurations, ' ...
    'all four branches = %.2e cm\n'], nTest, maxErr);
diary('off');
