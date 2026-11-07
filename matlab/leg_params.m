function P = leg_params(varargin)
%LEG_PARAMS Parameters of the 3-DoF robotic leg and of the stepping task.
%   P = LEG_PARAMS() returns the default values (lengths in cm, time in s).
%   P = LEG_PARAMS('R', 6, 'T', 3, ...) overrides any field by name.
%
%   Fields
%     l0, l1, l2, l3  link lengths (cm)
%     pA, pB          start and end of the step (cm), base frame {0}
%     T               period of one full step cycle (s); each of the two
%                     segments (semicircle, straight line) takes T/2
%     accFrac         fraction of each segment spent accelerating (and the
%                     same fraction decelerating) in the trapezoidal profile
%     dt              sampling period of the kinematic simulation (s)
%     cycles          number of step cycles to simulate
%
%   The semicircle has centre pC = (pA + pB)/2 and radius R = |pB - pA|/2,
%   and lies in the vertical plane through pA and pB (y0 points upwards).
%   The step runs along z0, so pA and pB may differ only in z.
%   'pC' and/or 'R' can be given instead of pA and pB.

P.l0 = 15;  P.l1 = 20;  P.l2 = 30;  P.l3 = 30;
P.pA = [20; -30; 45];
P.pB = [20; -30; 60];
P.T = 2;
P.accFrac = 1/3;
P.dt = 1e-3;
P.cycles = 1;

if mod(numel(varargin), 2) ~= 0
    error('leg_params:args', 'Parameters must be given as name-value pairs.');
end
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~isfield(P, name) && ~strcmp(name, 'R') && ~strcmp(name, 'pC')
        error('leg_params:unknown', 'Unknown parameter "%s".', name);
    end
    P.(name) = varargin{k+1};
end

% A centre and/or radius given by the user defines pA and pB along z0.
if isfield(P, 'pC') || isfield(P, 'R')
    if ~isfield(P, 'pC'), P.pC = (P.pA + P.pB)/2; end
    if ~isfield(P, 'R'),  P.R = norm(P.pB - P.pA)/2; end
    P.pC = P.pC(:);
    P.pA = P.pC - [0; 0; P.R];
    P.pB = P.pC + [0; 0; P.R];
end
P.pC = (P.pA + P.pB)/2;
P.R  = norm(P.pB - P.pA)/2;

assert(P.R > 0, 'The two stance points must differ.');
assert(abs(P.pA(2) - P.pB(2)) < 1e-12, 'pA and pB must lie on the same horizontal plane (equal y).');
assert(abs(P.pA(1) - P.pB(1)) < 1e-12, 'The step must run along z0: pA and pB must have equal x.');
assert(P.accFrac > 0 && P.accFrac <= 0.5, 'accFrac must be in (0, 0.5].');
end
