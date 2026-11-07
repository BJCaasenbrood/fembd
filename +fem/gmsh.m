function [status, out] = gmsh(args)
%GMSH Run the Gmsh executable with a command-line argument string.
%   [status, out] = fem.gmsh(args) uses $GMSH, otherwise 'gmsh' on the PATH
%   (on macOS also /Applications/Gmsh.app). Examples:
%     setenv('GMSH', 'C:\Program Files\gmsh-4.13.1-Windows64\gmsh.exe')
%     setenv('GMSH', '/usr/bin/gmsh')
%   Without output arguments the Gmsh output is printed, and an error is
%   thrown if Gmsh fails.

exe = strip(getenv('GMSH'), '"');
if isempty(exe)
    exe = 'gmsh';
    app = '/Applications/Gmsh.app/Contents/MacOS/gmsh';
    if ismac && isfile(app), exe = app; end
end
cmd = sprintf('"%s" %s', exe, args);
if isunix && ~ismac      % Linux: hide MATLAB's own (older) libstdc++ from Gmsh
    lib = strsplit(getenv('LD_LIBRARY_PATH'), ':');
    lib = lib(~startsWith(lib, matlabroot));
    cmd = sprintf('LD_LIBRARY_PATH="%s" %s', strjoin(lib, ':'), cmd);
end
[status, out] = system(cmd);
out = regexprep(out, '\x1b\[[0-9;]*m', '');      % strip colour codes
if nargout < 1
    if status ~= 0
        error('fem:gmsh', ['Gmsh failed:\n%s\nCheck that "%s" runs from a terminal, ' ...
              'or set the executable with setenv(''GMSH'', ''<path to gmsh>'').'], out, exe);
    end
    if ~isempty(strtrim(out)), disp(strtrim(out)); end
end
end
