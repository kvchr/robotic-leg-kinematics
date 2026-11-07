function [t, p, v, a, s] = trajectory_leg(P)
%TRAJECTORY_LEG Desired foot trajectory for one or more step cycles.
%   [t, p, v, a, s] = TRAJECTORY_LEG(P) returns the time vector t (1xN),
%   the desired foot position p, velocity v and acceleration a (3xN, cm,
%   cm/s, cm/s^2) and the path parameter s (1xN) of the current segment.
%
%   Each cycle of period T has two segments of duration T/2:
%     1) semicircle from pA over the top to pB (swing phase):
%          x = xA,  z = zC - R*cos(pi*s),  y = yC + R*sin(pi*s)
%     2) straight line from pB back to pA (stance phase):
%          x = xA,  z = zB + (zA - zB)*s,  y = yA
%   In both segments the path parameter s goes from 0 to 1 with a
%   trapezoidal speed profile: constant acceleration for accFrac*T/2,
%   constant speed, constant deceleration for the last accFrac*T/2.
%   The foot is at rest at pA and pB, so the velocity is continuous over
%   the whole cycle and from one cycle to the next.

Ts = P.T/2;
N = round(P.cycles*P.T/P.dt);
t = (0:N)*P.dt;
p = zeros(3, N+1); v = p; a = p; s = zeros(1, N+1);

for k = 1:N+1
    tc = mod(t(k), P.T);
    if k == N+1, tc = P.T; end              % end exactly at pA
    if tc <= Ts
        [sk, sd, sdd] = trapezoid(tc, Ts, P.accFrac);
        th = pi*sk;
        p(:,k) = [P.pA(1); P.pC(2) + P.R*sin(th); P.pC(3) - P.R*cos(th)];
        v(:,k) = pi*P.R*sd*[0; cos(th); sin(th)];
        a(:,k) = pi*P.R*sdd*[0; cos(th); sin(th)] + (pi*sd)^2*P.R*[0; -sin(th); cos(th)];
    else
        [sk, sd, sdd] = trapezoid(tc - Ts, Ts, P.accFrac);
        d = P.pA - P.pB;
        p(:,k) = P.pB + d*sk;
        v(:,k) = d*sd;
        a(:,k) = d*sdd;
    end
    s(k) = sk;
end
end

function [s, sd, sdd] = trapezoid(tau, Ts, f)
% Normalised trapezoidal profile on [0, Ts]: s(0) = 0, s(Ts) = 1,
% s'(0) = s'(Ts) = 0, with acceleration time ta = f*Ts.
ta = f*Ts;
vmax = 1/(Ts - ta);
acc = vmax/ta;
if tau < ta
    s = 0.5*acc*tau^2;            sd = acc*tau;     sdd = acc;
elseif tau <= Ts - ta
    s = 0.5*acc*ta^2 + vmax*(tau - ta);  sd = vmax;  sdd = 0;
else
    r = Ts - tau;
    s = 1 - 0.5*acc*r^2;          sd = acc*r;       sdd = -acc;
end
end
