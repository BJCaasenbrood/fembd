# fmbd — Flexible Multi-Body Dynamics in MATLAB

fmbd builds dynamic models of mechanisms from STL parts. Rigid parts only need their surface or
their mass properties. Elastic parts are meshed with Gmsh and reduced to a few modes that move with
the part. The parts are connected with joints, couplings, ties, springs, dampers and
transmissions such as belts and gears.

A model can be simulated through large motions with a controller in the loop, solved for its
static equilibrium, and linearised about any pose to get modes and state-space models. The linear
models include the stiffness from gravity and preloads, so they change with the pose of the
mechanism, as they do for a real robot arm.

```matlab
addpath('path/to/PDEToolbox');
setenv('GMSH', '/path/to/gmsh'); % if gmsh is not on the PATH

alu  = fmbd.material(70e9, 0.33, 2700);
link = fmbd.Body.rigid('link.stl', 7850, Scale = 1e-3); % STL in mm
tool = fmbd.Body.elastic('tool.stl', alu, Size = 3, Scale = 1e-3, Modes = 10);

s = fmbd.System(Gravity = [0 0 -9.81]);
s.add(link);
s.add(tool);
s.joint('ground', 'link', 'revolute', Point = [0 0 0], Name = 'shoulder', Position = 0.3);
s.joint('link', s.couple('tool', @(X) X(:,1) < 0.201), 'fixed'); % bolt the tool's flange
s.input('tau', 'shoulder');
s.output('th', 'shoulder');
s.output('w', 'shoulder', Type = 'velocity');

% a PD-controlled move
[t, y, x] = s.simulate([0 1], Input = @(t, y) 80*(1 - y.th) - 8*y.w);
s.animate(t, x, Scale = 100); % deformation magnified

% linear analysis about the equilibrium
x0 = s.equilibrium();
[f, V] = s.modes(x0); % s.plotMode(x0, V(:,k)) shows mode k
sys = s.statespace(x0, Modes = 20, Format = "ss"); % needs the Control System Toolbox
```

## Installing Gmsh
fmbd calls the Gmsh executable (version 4). It uses `$GMSH` if that is set, and otherwise `gmsh`
on the PATH.

- Windows: download the 64-bit zip from https://gmsh.info, unzip it, and set
  `setenv('GMSH', 'C:\path\to\gmsh-4.x.x-Windows64\gmsh.exe')`.
- Linux: `sudo apt install gmsh`, or the binary from gmsh.info. MATLAB's own `libstdc++` is kept
  away from Gmsh automatically.

Put the `setenv` line in your `startup.m`, and check the setup with `fmbd.gmsh('--version')`.

## Getting started
The [`tutorials`](tutorials/README.md) folder has eight scripts that build up from a first model
to elastic bodies in motion. The [`examples`](examples) folder has three complete models: an STL
assembly with a rigid payload, a two-link arm at several poses, and a parallel leaf-spring guide.
All functions have help text, e.g. `help fmbd.System/simulate`.

fmbd needs MATLAB R2022a or later and Gmsh 4.
