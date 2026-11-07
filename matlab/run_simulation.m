function S = run_simulation(varargin)
%RUN_SIMULATION Kinematic simulation of the 3-DoF leg (Part B).
%   S = RUN_SIMULATION() runs the default task and saves the figures in
%   ../figures and a summary in ../results/simulation_summary.txt.
%   S = RUN_SIMULATION('R', 6, 'T', 3, ...) changes task parameters
%   (see leg_params). Every sample is checked: the run stops with an
%   error if a point is unreachable or the leg comes close to a singular
%   configuration (normalized det J below minDet = 0.01).
%   RUN_SIMULATION(..., 'save', false) only shows the figures.
%
%   Steps: desired foot trajectory -> inverse kinematics q(t) ->
%   joint velocities qdot = J^-1 * v -> checks (forward kinematics,
%   numerical derivative, distance from singular configurations) ->
%   plots, stick diagrams and an animated GIF.

doSave = true;
k = find(strcmp(varargin, 'save'), 1);
if ~isempty(k), doSave = varargin{k+1}; varargin(k:k+1) = []; end
P = leg_params(varargin{:});

here = fileparts(mfilename('fullpath'));
figDir = fullfile(here, '..', 'figures');
resDir = fullfile(here, '..', 'results');
if doSave
    if ~exist(figDir, 'dir'), mkdir(figDir); end
    if ~exist(resDir, 'dir'), mkdir(resDir); end
end

% ---------- 1. Desired trajectory in task space ----------
[t, p, v, a] = trajectory_leg(P);
N = numel(t);

% ---------- 2. Inverse kinematics and joint velocities ----------
minDet = 0.01;                              % smallest allowed normalized |det J|
detScale = P.l2*P.l3*(P.l2 + P.l3);
q = zeros(3, N); qd = q; detJ = zeros(1, N); fkErr = zeros(1, N);
for k = 1:N
    q(:,k) = ik_leg(p(:,k), P).';
    [JL, detJ(k)] = jacobian_leg(q(:,k), P);
    if abs(detJ(k))/detScale < minDet
        error('run_simulation:singular', ['The leg is too close to a singular configuration at ' ...
            't = %.3f s (normalized det J = %.4f).'], t(k), detJ(k)/detScale);
    end
    qd(:,k) = JL \ v(:,k);
    fkErr(k) = norm(fk_leg(q(:,k), P) - p(:,k));
end
q = unwrap(q, [], 2);                       % continuous angles
qdNum = gradient(q, P.dt);                  % numerical check of J^-1 v
qdErr = max(abs(qdNum(:, 2:end-1) - qd(:, 2:end-1)), [], 'all');
detRel = detJ/detScale;                     % dimensionless distance from singularity

% ---------- 3. Summary ----------
speed = vecnorm(v);
Ts = P.T/2;
swing = mod(t, P.T) <= Ts;                  % samples of the semicircle segments
qdeg = rad2deg(q);  qdeg(abs(qdeg) < 5e-7) = 0;   % print 0 instead of -0
qdmax = rad2deg(max(abs(qd), [], 2));  qdmax(qdmax < 5e-7) = 0;
S = struct('t', t, 'p', p, 'v', v, 'a', a, 'q', q, 'qd', qd, 'detJ', detJ, 'P', P);
lines = {
    sprintf('Robotic leg kinematic simulation (MATLAB %s)', version('-release'))
    sprintf('Links l0 = %g, l1 = %g, l2 = %g, l3 = %g cm', P.l0, P.l1, P.l2, P.l3)
    sprintf('pA = (%g, %g, %g) cm, pB = (%g, %g, %g) cm, centre pC = (%g, %g, %g) cm, R = %g cm', P.pA, P.pB, P.pC, P.R)
    sprintf('Period T = %g s (two segments of %g s), acceleration fraction %.3f, dt = %g s, cycles = %d', P.T, Ts, P.accFrac, P.dt, P.cycles)
    sprintf('Path length: semicircle %.3f cm, line %.3f cm', pi*P.R, 2*P.R)
    sprintf('Peak foot speed: swing %.3f cm/s, stance %.3f cm/s (theory %.3f and %.3f)', ...
        max(speed(swing)), max(speed(~swing)), pi*P.R/(Ts*(1-P.accFrac)), 2*P.R/(Ts*(1-P.accFrac)))
    sprintf('Foot height above the walking plane: max %.3f cm', max(p(2,:)) - P.pA(2))
    sprintf('Max |FK(q(t)) - p_d(t)| = %.2e cm', max(fkErr))
    sprintf('Max |J^-1 v - dq/dt (numerical)| = %.2e rad/s', qdErr)
    sprintf('q1 range [%.3f, %.3f] deg', min(qdeg(1,:)), max(qdeg(1,:)))
    sprintf('q2 range [%.3f, %.3f] deg', min(qdeg(2,:)), max(qdeg(2,:)))
    sprintf('q3 range [%.3f, %.3f] deg', min(qdeg(3,:)), max(qdeg(3,:)))
    sprintf('Peak |qdot|: q1 %.3f, q2 %.3f, q3 %.3f deg/s', qdmax)
    sprintf('Min |cos q3| = %.3f, min |l3*c23 - l2*s2| = %.3f cm, min |det J|/(l2*l3*(l2+l3)) = %.3f', ...
        min(abs(cos(q(3,:)))), min(abs(detJ./(P.l2*P.l3*cos(q(3,:))))), min(abs(detRel)))
    };
fprintf('%s\n', lines{:});
if doSave
    fid = fopen(fullfile(resDir, 'simulation_summary.txt'), 'w');
    fprintf(fid, '%s\n', lines{:});
    fclose(fid);
end

% ---------- 4. Figures ----------
f0 = plot_frames(P);
save_fig(f0, figDir, 'fig0_dh_frames.png', doSave);

lab = {'x', 'y', 'z'};
cols = lines_colors();

f1 = newfig([900 520]);
tl = tiledlayout(f1, 2, 1, 'TileSpacing', 'compact');
nexttile(tl); hold on; grid on;
for i = 1:3, plot(t, p(i,:), 'LineWidth', 1.6, 'Color', cols(i,:)); end
xline(Ts*(1:2*P.cycles-1), ':'); ylabel('position (cm)');
legend(strcat('p_{E', lab, '}'), 'Location', 'eastoutside');
title('Desired foot position');
nexttile(tl); hold on; grid on;
for i = 1:3, plot(t, v(i,:), 'LineWidth', 1.6, 'Color', cols(i,:)); end
plot(t, speed, 'k--', 'LineWidth', 1.2);
xline(Ts*(1:2*P.cycles-1), ':'); ylabel('velocity (cm/s)'); xlabel('time t (s)');
legend([strcat('v_{E', lab, '}'), {'|v_E|'}], 'Location', 'eastoutside');
title('Desired foot velocity');
save_fig(f1, figDir, 'fig1_task_space_profile.png', doSave);

f2 = newfig([900 520]);
tl = tiledlayout(f2, 2, 1, 'TileSpacing', 'compact');
nexttile(tl); hold on; grid on;
for i = 1:3, plot(t, rad2deg(q(i,:)), 'LineWidth', 1.6, 'Color', cols(i,:)); end
xline(Ts*(1:2*P.cycles-1), ':'); ylabel('angle (deg)');
legend({'$q_1$', '$q_2$', '$q_3$'}, 'Interpreter', 'latex', 'FontSize', 11, 'Location', 'eastoutside');
title('Joint angles');
nexttile(tl); hold on; grid on;
for i = 1:3, plot(t, rad2deg(qd(i,:)), 'LineWidth', 1.6, 'Color', cols(i,:)); end
idx = 1:40:N;
for i = 1:3, plot(t(idx), rad2deg(qdNum(i,idx)), 'o', 'MarkerSize', 3, 'Color', cols(i,:)); end
xline(Ts*(1:2*P.cycles-1), ':'); ylabel('angular velocity (deg/s)'); xlabel('time t (s)');
legend({'$\dot{q}_1$', '$\dot{q}_2$', '$\dot{q}_3$', 'numerical $dq/dt$'}, 'Interpreter', 'latex', 'FontSize', 11, 'Location', 'eastoutside');
title('Joint velocities: J^{-1} v (lines) and numerical derivative (markers)');
save_fig(f2, figDir, 'fig2_joint_space_profile.png', doSave);

% Stick diagrams: 3D view and side view (z0 horizontal, y0 vertical).
nSnap = 11;
snapIdx = round(linspace(1, round(P.T/P.dt) + 1, nSnap));
f3 = newfig([1000 460]);
tl = tiledlayout(f3, 1, 2, 'TileSpacing', 'compact');
ax1 = nexttile(tl); hold(ax1, 'on'); grid(ax1, 'on');
ax2 = nexttile(tl); hold(ax2, 'on'); grid(ax2, 'on');
cm = parula(nSnap + 2);
for j = 1:nSnap
    [~, J] = fk_leg(q(:, snapIdx(j)), P);
    plot3(ax1, J(1,:), J(3,:), J(2,:), '-o', 'Color', cm(j,:), 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', cm(j,:));
    plot(ax2, J(3,:), J(2,:), '-o', 'Color', cm(j,:), 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', cm(j,:));
end
plot3(ax1, p(1,:), p(3,:), p(2,:), 'r-', 'LineWidth', 2);
plot(ax2, p(3,:), p(2,:), 'r-', 'LineWidth', 2);
plot3(ax1, 0, 0, 0, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
plot(ax2, 0, 0, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
xlabel(ax1, 'x_0 (cm)'); ylabel(ax1, 'z_0 (cm)'); zlabel(ax1, 'y_0 (cm)');
axis(ax1, 'equal'); axis(ax1, [-5 30 -5 65 -35 10]); view(ax1, -40, 18);
title(ax1, '3D view (y_0 points up)');
xlabel(ax2, 'z_0 (cm)'); ylabel(ax2, 'y_0 (cm)'); axis(ax2, 'equal');
title(ax2, 'Side view (plane of the step)');
cb = colorbar(ax2); colormap(f3, cm(1:nSnap,:)); clim(ax2, [0 P.T]); cb.Label.String = 'time t (s)';
title(tl, sprintf('Leg configurations every %.1f s over one cycle (red: foot path)', P.T/(nSnap-1)));
save_fig(f3, figDir, 'fig3_stick_diagram.png', doSave);

f4 = newfig([900 300]);
plot(t, detRel, 'LineWidth', 1.6); grid on; hold on;
yline(0, 'k-'); xline(Ts*(1:2*P.cycles-1), ':');
xlabel('time t (s)'); ylabel('det J_L / (l_2 l_3 (l_2+l_3))');
title('Distance from singular configurations (zero means singular)');
save_fig(f4, figDir, 'fig4_jacobian_determinant.png', doSave);

% Animated GIF of one cycle (side and 3D view).
if doSave
    gifFile = fullfile(figDir, 'leg_animation.gif');
    if exist(gifFile, 'file'), delete(gifFile); end
    f5 = newfig([720 360]);
    tl = tiledlayout(f5, 1, 2, 'TileSpacing', 'compact');
    b1 = nexttile(tl); b2 = nexttile(tl);
    frames = 1:40:round(P.T/P.dt) + 1;
    for j = frames
        [~, J] = fk_leg(q(:, j), P);
        cla(b1); hold(b1, 'on'); grid(b1, 'on');
        plot3(b1, p(1,:), p(3,:), p(2,:), 'r-', 'LineWidth', 1.2);
        plot3(b1, J(1,:), J(3,:), J(2,:), '-o', 'Color', [0 0.45 0.74], 'LineWidth', 2.5, 'MarkerFaceColor', [0 0.45 0.74]);
        axis(b1, 'equal'); axis(b1, [-5 30 -5 65 -35 10]); view(b1, -40, 18);
        xlabel(b1, 'x_0'); ylabel(b1, 'z_0'); zlabel(b1, 'y_0');
        cla(b2); hold(b2, 'on'); grid(b2, 'on');
        plot(b2, p(3,:), p(2,:), 'r-', 'LineWidth', 1.2);
        plot(b2, J(3,:), J(2,:), '-o', 'Color', [0 0.45 0.74], 'LineWidth', 2.5, 'MarkerFaceColor', [0 0.45 0.74]);
        axis(b2, 'equal'); axis(b2, [-5 65 -35 10]); xlabel(b2, 'z_0'); ylabel(b2, 'y_0');
        title(tl, sprintf('t = %.2f s', t(j)));
        plot3(b1, 0, 0, 0, 'ks', 'MarkerFaceColor', 'k');
        plot(b2, 0, 0, 'ks', 'MarkerFaceColor', 'k');
        exportgraphics(f5, gifFile, 'Append', true, 'Resolution', 110);
    end
    close(f5);
end
end

% ---------- helpers ----------
function f = newfig(sz)
f = figure('Color', 'w', 'Position', [100 100 sz]);
end

function save_fig(f, d, name, doSave)
if doSave
    exportgraphics(f, fullfile(d, name), 'Resolution', 150);
    close(f);
end
end

function c = lines_colors()
c = [0 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188];
end
