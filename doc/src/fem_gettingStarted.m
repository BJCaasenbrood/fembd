%% Getting Started
%
%% Installation
% Add the toolbox folder (the folder containing the |+fem| package and
% |info.xml|) to the MATLAB path, for example in |startup.m|:
%
%   addpath('C:\path\to\PDEToolbox')
%
% Requirements: MATLAB R2022a or later, and Gmsh 4 for meshing STL files.
%
%% Gmsh
% Bodies are meshed by calling the Gmsh executable. fem uses the GMSH
% environment variable if it is set, otherwise |gmsh| on the system path.
%
% * *Windows*: download the Windows 64-bit zip from
%   <https://gmsh.info gmsh.info>, unzip it anywhere and set
%   |setenv('GMSH', 'C:\gmsh-4.x.x-Windows64\gmsh.exe')|.
% * *Linux*: |sudo apt install gmsh|. MATLAB's own |libstdc++| is
%   hidden from Gmsh automatically.
% * *macOS*: install Gmsh.app in /Applications (found automatically).
%
% Put the |setenv| line in |startup.m|. Check the setup with:
%
%   fem.gmsh('--version')
%
% Rigid bodies created with |fem.Body.rigidSurface| or |fem.Body.pointMass|
% do not need Gmsh.
%
%% A first model
% A 1000 x 50 x 50 mm steel cantilever. The STL is in millimetres, so
% |Scale = 1e-3| converts the mesh to metres; |Size| stays in file units.
fem.boxSTL('beam.stl', [0 0 0], [1000 50 50]);
steel = fem.material(210e9, 0.3, 7850);              % E, nu, rho
beam  = fem.Body('beam.stl', steel, Size = 25, Scale = 1e-3);

a = fem.Assembly();
a.add(beam);
a.fix('beam', @(X) X(:,1) < 1e-9);                    % clamp x = 0

[K, M, T] = a.system();                               % reduced K, M; u = T*q
[V, f] = fem.modes(K, M, 6);                          % Hz, mass-normalised
a.plot(T*V(:,1), 0.1);                                % mode 1, scaled
%%
% The workflow has three stages:
%
% # *Bodies* (|fem.Body|): mesh, material, and the body's own K and M.
% # *Assembly* (|fem.Assembly|): bodies plus connections and supports.
% # *System* (|a.system|): constraints are eliminated, u = T*q, giving
%   the reduced matrices K = T'*Kfull*T and M = T'*Mfull*T.
%
%% Next steps
%
% * <tut01_firstModel.html Tutorial 1> walks through this model and checks
%   it against beam theory.
% * <fem_tutorials.html All tutorials>
% * <fem_reference.html Function reference>
