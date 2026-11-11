function setVertices(h, V)
%SETVERTICES Move bodies drawn by drawBodies to the positions V{i}.
for i = 1:numel(h)
    if isa(h(i), 'matlab.graphics.primitive.Patch')
        set(h(i), 'Vertices', V{i});
    else
        set(h(i), 'XData', V{i}(1), 'YData', V{i}(2), 'ZData', V{i}(3));
    end
end
end
