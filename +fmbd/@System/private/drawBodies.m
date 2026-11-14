function h = drawBodies(bodies, V, C)
% Draw bodies as patches (point masses as markers) at node positions V{i}
% with colours C{i}; returns one handle per body.
h = gobjects(1, numel(bodies));
held = ishold;  hold on
for i = 1:numel(bodies)
    F = bodies{i}.F;
    if isempty(F)
        h(i) = plot3(V{i}(1), V{i}(2), V{i}(3), 'ko', 'MarkerFaceColor', C{i}(1,:), 'MarkerSize', 8);
        continue
    end
    if size(F, 2) == 6 % Tet10 faces: split into 4 triangles
        F = [F(:,[1 4 6]); F(:,[4 2 5]); F(:,[6 5 3]); F(:,[4 5 6])];
    end
    Ci = C{i};  if size(Ci, 1) == 1, Ci = repmat(Ci, size(V{i}, 1), 1); end
    h(i) = patch('Faces', F, 'Vertices', V{i}, 'FaceVertexCData', Ci, ...
                 'FaceColor', 'interp', 'LineStyle', 'none');
end
if ~held, hold off, end
axis equal;
view(3);
shading interp
if isempty(findobj(gca, 'Type', 'light')), camlight, end
end
