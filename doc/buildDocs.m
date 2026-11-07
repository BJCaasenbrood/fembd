function buildDocs(opts)
%BUILDDOCS Build the fem documentation for the MATLAB Help browser.
%   buildDocs publishes the pages in doc/src, the tutorials and the examples
%   to doc/html, copies the table of contents, and builds the search index.
%   Tutorials and examples are run (Gmsh is needed) so that their output and
%   figures appear in the pages; this takes a few minutes.
%
%   Options: Run (default true) - run tutorials and examples while
%            publishing; false publishes their code only.
%
%   After building, the documentation appears in the Help browser under
%   Supplemental Software (the toolbox folder, which contains info.xml,
%   must be on the path).

arguments
    opts.Run (1,1) logical = true
end
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
out  = fullfile(here, 'html');
addpath(root);
if ~isfolder(out), mkdir(out); end

% reference pages: code is shown, not run
for f = dir(fullfile(here, 'src', '*.m'))'
    publish(fullfile(f.folder, f.name), outputDir = out, evalCode = false);
end

% tutorials and examples, with output and figures
scripts = [dir(fullfile(root, 'tutorials', 'tut*.m')); dir(fullfile(root, 'examples', 'demo*.m'))]';
for f = scripts
    fprintf('publishing %s\n', f.name);
    evalIt = opts.Run && ~strcmp(f.name, 'demoParallelFlexure.m');   % code only
    publishScript(fullfile(f.folder, f.name), out, evalIt);
end

renderEquations(out);
copyfile(fullfile(here, 'helptoc.xml'), out);
builddocsearchdb(out);
fprintf('Documentation built in %s\n', out);
end

function publishScript(file, out, run)
% publish evaluates code in the base workspace: start each script clean.
evalin('base', 'clear');  close all
folder = fileparts(file);  addpath(folder);            % publish needs it on the path
cleanup = onCleanup(@() rmpath(folder));
publish(file, outputDir = out, evalCode = run, catchError = false);
close all
end

function renderEquations(out)
% publish renders equations at screen resolution (tiny without a display);
% re-render them from their LaTeX at 150 dpi and size them accordingly.
dpi = 150;
for f = dir(fullfile(out, '*.html'))'
    file = fullfile(out, f.name);  html = fileread(file);
    tok = regexp(html, '<img src="([^"]+_eq\d+\.png)" alt="([^"]*)"[^>]*>', 'tokens');
    for t = tok
        [img, tex] = t{1}{:};
        tex = replace(tex, {'&#10;', '&lt;', '&gt;', '&quot;', '&amp;'}, {' ', '<', '>', '"', '&'});
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 1200 300]);
        ax = axes(fig, 'Position', [0 0 1 1], 'Visible', 'off');
        text(ax, 0.5, 0.5, tex, 'Interpreter', 'latex', 'FontSize', 13, ...
             'HorizontalAlignment', 'center');
        exportgraphics(fig, fullfile(out, img), 'Resolution', dpi);
        close(fig);
        sz = size(imread(fullfile(out, img)));  css = round(sz([2 1])*96/dpi);
        html = regexprep(html, ['(<img src="' regexptranslate('escape', img) '"[^>]*?)\s*style="[^"]*"'], ...
                         sprintf('$1 style="width:%dpx;height:%dpx;"', css));
    end
    writelines(html, file);
end
end
