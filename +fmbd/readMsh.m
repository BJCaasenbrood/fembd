function [X, T] = readMsh(file)
%READMSH Read tetrahedra from an ASCII Gmsh mesh file (format 2.2 or 4.1).
%   [X, T] = fmbd.readMsh(file) returns nodes X (n x 3) and tetrahedra T
%   (ne x 4 or ne x 10, Gmsh node ordering). Only the highest-order tet
%   type present is returned; unused nodes are removed.

txt = string(fileread(file));
sec = @(name) sscanf(extractBetween(txt, "$" + name, "$End" + name), '%f');
fmt = sec('MeshFormat');
assert(fmt(2) == 0, 'fmbd:readMsh', 'Binary .msh files are not supported.');

if fmt(1) < 3                                   % ---- MSH 2.2
    v = sec('Nodes');
    v = reshape(v(2:end), 4, []);
    tag = v(1,:);  xyz = v(2:4,:)';
    lines = cellstr(splitlines(strtrim(extractBetween(txt, "$Elements", "$EndElements"))));
    e = cellfun(@(l) sscanf(l, '%f')', lines(2:end), 'UniformOutput', false);
    type = cellfun(@(r) r(2), e);
    tet = max(type(type == 4 | type == 11));
    assert(~isempty(tet), 'fmbd:readMsh', 'No tetrahedra in %s.', file);
    e = vertcat(e{type == tet});
    T = e(:, end - (tet == 11)*6 - 3 : end);
else                                            % ---- MSH 4.1
    v = sec('Nodes');
    xyz = zeros(v(2), 3);  tag = zeros(v(2), 1);
    p = 5;  k = 0;
    for b = 1:v(1)
        [dim, param, n] = deal(v(p), v(p+2), v(p+3));  p = p + 4;
        nc = 3 + param*dim;
        tag(k+1:k+n) = v(p:p+n-1);  p = p + n;
        c = reshape(v(p:p+nc*n-1), nc, n)';  p = p + nc*n;
        xyz(k+1:k+n, :) = c(:, 1:3);  k = k + n;
    end
    e = sec('Elements');
    nper = [2 3 4 4 8 6 5 3 6 9 10 27 18 14 1 8 20 15 13];  % nodes per type
    blk = {};  p = 5;
    for b = 1:e(1)
        [type, n] = deal(e(p+2), e(p+3));  p = p + 4;
        m = nper(type) + 1;
        if type == 4 || type == 11
            blk(end+1, :) = {type, reshape(e(p:p+m*n-1), m, n)'}; %#ok<AGROW>
        end
        p = p + m*n;
    end
    assert(~isempty(blk), 'fmbd:readMsh', 'No tetrahedra in %s.', file);
    tet = max([blk{:,1}]);
    T = vertcat(blk{[blk{:,1}] == tet, 2});
    T = T(:, 2:end);
end

map(tag) = 1:numel(tag);
T = map(T);
[used, ~, ic] = unique(T);
X = xyz(used, :);
T = reshape(ic, size(T));
end
