# fembd — finite elements and multibody dynamics in MATLAB

Rigid and elastic multibody systems from STL assemblies. Parts are meshed with Gmsh, given
linear materials or mass properties, and connected by joints, couplings, ties, springs, dampers
and transmissions. One model gives two analyses:

- **Nonlinear simulation** (`a.nonlinear`): large motions with exact joint kinematics, gravity,
  controllers in the loop, and elastic bodies in a floating frame (Craig–Bampton modes).
- **Linear analysis about a pose** (`a.system`, `a.statespace`): the full finite-element model
  for modes, frequency responses, statics, state space for control, and model reduction.

```matlab
addpath('path/to/PDEToolbox');
setenv('GMSH', '/path/to/gmsh');                 % if gmsh is not on the PATH

alu  = fembd.material(70e9, 0.33, 2700);
link = fembd.Body.rigidSurface('link.stl', 7850, Scale = 1e-3);    % rigid, STL in mm
tool = fembd.Body('tool.stl', alu, Size = 3, Scale = 1e-3);        % elastic, Gmsh mesh

a = fembd.Assembly();
a.add(link);  a.add(tool);
a.joint('ground', 'link', 'revolute', Point = [0 0 0], Name = 'shoulder');
a.joint('link', a.couple('tool', @(X) X(:,1) < 0.201), 'fixed');  % bolt the tool's flange
a.input('tau', 'shoulder');
a.output('th', 'shoulder');  a.output('w', 'shoulder', Type = 'velocity');

% nonlinear: a PD-controlled move under gravity
m = a.nonlinear(Gravity = [0 0 -9.81]);
[t, y, x] = m.simulate([0 1], Input = @(t, y) 80*(1 - y.th) - 8*y.w);
m.animate(t, x, Scale = 100);                    % deformation magnified

% linear: modes and a state-space model about the current pose
[K, M, T] = a.system();                          % reduced matrices, u = T*q
[V, f] = fembd.modes(K, M, 10);                  % Hz, mass-normalised
sys = a.statespace(Modes = 20, Format = "ss");   % tau -> th, w (Control System Toolbox)
```

## Installing Gmsh
The toolbox calls the Gmsh executable (version 4). It uses `$GMSH` if that is set, otherwise
`gmsh` on the PATH.
- **Windows:** download the Windows 64-bit zip from https://gmsh.info, unzip it anywhere, then
  `setenv('GMSH', 'C:\path\to\gmsh-4.x.x-Windows64\gmsh.exe')`.
- **Linux:** `sudo apt install gmsh` (or the gmsh.info binary). MATLAB's own `libstdc++` is hidden
  from Gmsh automatically.
- **macOS:** install `Gmsh.app` in `/Applications` (found automatically), or `brew install gmsh`.

Put the `setenv` line in `startup.m` to keep it across sessions. Run `fembd.gmsh('--version')` to
check the setup.

## Documentation
- **Tutorials:** [`tutorials/`](tutorials/README.md), ten runnable scripts: a first nonlinear model, then meshing, materials, rigid bodies, joints and ties, then linear analysis (modes, statics, dynamics, reduction).
- **Manual:** [`docs/fembd.md`](docs/fembd.md): workflow, the nonlinear model, the formulation, the API, linear-analysis recipes and limitations.
- **In MATLAB:** `doc fembd.Assembly`, `doc fembd.Model`, `help fembd`. The tutorials, examples and function reference also appear in the Help browser under *Supplemental Software* (prebuilt in `doc/html`; rebuild with `cd doc; buildDocs`, which runs the tutorials and needs Gmsh).
- **Examples:** [`examples/`](examples): an STL assembly with a rigid payload, a two-link arm on joints, and a parallel leaf-spring guide.
- **Verification:** `tests/runTests.m`.

Requires MATLAB R2022a or later and Gmsh 4.
