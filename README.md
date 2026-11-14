# fmbd — flexible multibody dynamics in MATLAB

Rigid and elastic multibody systems from STL assemblies. Rigid parts need only their surface or
mass properties; elastic parts are meshed with Gmsh and reduced to a few modes in a moving frame.
Parts are connected by joints, couplings, ties, springs, dampers and transmissions. One nonlinear
model gives:

- **simulation** of large motions with exact joint kinematics, gravity and controllers in the loop;
- **equilibria** under gravity, springs and constant inputs;
- **linear models** about any operating point (modes, state space, frequency responses), with the
  stiffness from gravity, constraint forces and preload included.

```matlab
addpath('path/to/PDEToolbox');
setenv('GMSH', '/path/to/gmsh');                 % if gmsh is not on the PATH

alu  = fmbd.material(70e9, 0.33, 2700);
link = fmbd.Body.rigid('link.stl', 7850, Scale = 1e-3);                     % rigid, STL in mm
tool = fmbd.Body.elastic('tool.stl', alu, Size = 3, Scale = 1e-3, Modes = 10);  % elastic

s = fmbd.System(Gravity = [0 0 -9.81]);
s.add(link);  s.add(tool);
s.joint('ground', 'link', 'revolute', Point = [0 0 0], Name = 'shoulder', Position = 0.3);
s.joint('link', s.couple('tool', @(X) X(:,1) < 0.201), 'fixed');   % bolt the tool's flange
s.input('tau', 'shoulder');
s.output('th', 'shoulder');  s.output('w', 'shoulder', Type = 'velocity');

% simulation: a PD-controlled move
[t, y, x] = s.simulate([0 1], Input = @(t, y) 80*(1 - y.th) - 8*y.w);
s.animate(t, x, Scale = 100);                    % deformation magnified

% linear analysis about the equilibrium
x0 = s.equilibrium();
[f, V] = s.modes(x0);                            % Hz;  s.plotMode(x0, V(:,k))
sys = s.statespace(x0, Modes = 20, Format = "ss");   % tau -> th, w (Control System Toolbox)
```

## Installing Gmsh
The toolbox calls the Gmsh executable (version 4). It uses `$GMSH` if that is set, otherwise
`gmsh` on the PATH.
- **Windows:** download the Windows 64-bit zip from https://gmsh.info, unzip it anywhere, then
  `setenv('GMSH', 'C:\path\to\gmsh-4.x.x-Windows64\gmsh.exe')`.
- **Linux:** `sudo apt install gmsh` (or the gmsh.info binary). MATLAB's own `libstdc++` is hidden
  from Gmsh automatically.
- **macOS:** install `Gmsh.app` in `/Applications` (found automatically), or `brew install gmsh`.

Put the `setenv` line in `startup.m` to keep it across sessions. Run `fmbd.gmsh('--version')` to
check the setup.

## Documentation
- **Tutorials:** [`tutorials/`](tutorials/README.md), eight runnable scripts: a first model, meshing, materials, rigid bodies, joints and transmissions, ties and supports, linear analysis, elastic bodies in motion.
- **In MATLAB:** `help fmbd`, `help fmbd.System`, `help fmbd.System/simulate`, `doc fmbd.System`, etc.
- **Examples:** [`examples/`](examples): an STL assembly with a rigid payload, a two-link arm at several poses, and a parallel leaf-spring guide.

Requires MATLAB R2022a or later and Gmsh 4.
