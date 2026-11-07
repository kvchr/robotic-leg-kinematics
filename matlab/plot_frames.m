function f = plot_frames(P, q)
%PLOT_FRAMES Draw the leg with its Denavit-Hartenberg frames.
%   f = PLOT_FRAMES(P, q) draws the leg at joint angles q (default zeros)
%   and the axes x (red), y (green), z (blue) of the frames {0}, {1}, {2},
%   {3}, {E}. The plot uses x0 and z0 as horizontal axes and y0 upwards.

if nargin < 2, q = [0 0 0]; end
[~, J, ~, F] = fk_leg(q, P);
names = {'0', '1', '2', '3', 'E'};
L = 7;                                     % axis arrow length (cm)
colors = [0.85 0.1 0.1; 0.1 0.6 0.1; 0.1 0.2 0.85];
axl = 'xyz';

f = figure('Color', 'w', 'Position', [100 100 760 560]);
ax = axes(f); hold(ax, 'on'); grid(ax, 'on'); ax.Toolbar.Visible = 'off';
plot3(ax, J(1,:), J(3,:), J(2,:), '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 5);
plot3(ax, J(1,:), J(3,:), J(2,:), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
% Small offsets keep coincident frames ({2} and {3}) readable.
shift = {[0;0;0], [0;0;0], [0;0;0], [6;0;0], [0;0;0]};
% At q = 0, z0 and y1 both lie along the link l0: move their labels to
% opposite sides of it (offsets in base coordinates [x; y; z], y up).
nudge = containers.Map({'z_0', 'y_1'}, {[0; -2.5; 4], [0; 3; 0]});
for i = 1:numel(F)
    o = F{i}(1:3,4) + shift{i};
    for k = 1:3
        d = L*F{i}(1:3,k);
        quiver3(ax, o(1), o(3), o(2), d(1), d(3), d(2), 0, 'Color', colors(k,:), ...
            'LineWidth', 1.8, 'MaxHeadSize', 0.6);
        name = sprintf('%s_%s', axl(k), names{i});
        pos = o + 1.15*d;
        if isKey(nudge, name), pos = pos + nudge(name); end
        text(ax, pos(1), pos(3), pos(2), name, 'Color', colors(k,:), 'FontWeight', 'bold');
    end
    text(ax, o(1)+1, o(3)+1, o(2)-4, sprintf('O_%s', names{i}), 'FontSize', 11);
end
labels = {'l_0', 'l_1', 'l_2', 'l_3'};
lofs = {[2; -7; 0], [1; 3; 1], [1; 3; 1], [1; 3; 1]};   % l0 label below its link
for k = 1:4
    m = (J(:,k) + J(:,k+1))/2 + lofs{k};
    text(ax, m(1), m(3), m(2), labels{k}, 'FontSize', 12, 'FontAngle', 'italic');
end
axis(ax, 'equal'); view(ax, 55, 20);
xlabel(ax, 'x_0 (cm)'); ylabel(ax, 'z_0 (cm)'); zlabel(ax, 'y_0 (cm)');
title(ax, sprintf('D-H frames at q = (%g, %g, %g) deg', rad2deg(q)));
end
