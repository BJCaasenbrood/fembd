# fembd — finite elements and multibody dynamics in MATLAB

A compact research alternative to the PDE Toolbox for linear elasticity. It loads STL parts,
meshes them with Gmsh, assigns linear materials, connects the bodies, and returns the stiffness
and mass matrices for dynamic analysis. The same assembly also gives a state-space model for
control, and a nonlinear multibody model for large rotations (`a.nonlinear`).

```matlab
addpath('path/to/PDEToolbox');
setenv('GMSH', '/path/to/gmsh');                 % if gmsh is not on the PATH

steel = fembd.material(210e9, 0.30, 7850);
alu   = fembd.material(70e9,  0.33, 2700);
base   = fembd.Body('base.stl',   steel, Size = 8, Scale = 1e-3);   % STL in mm
pillar = fembd.Body('pillar.stl', alu,   Size = 3, Scale = 1e-3);

a = fembd.Assembly();
a.add(base);  a.add(pillar);
a.tie('base', 'pillar');                         % bonded interface (master, slave)
a.fix('base', @(X) X(:,3) < 1e-9);               % clamp nodes selected by position

[K, M, T] = a.system();                          % reduced matrices, u = T*q
[V, f] = fembd.modes(K, M, 10);                    % Hz, mass-normalised
a.plot(T*V(:,1), 0.01);

ib = a.qdof(a.dof('pillar', @(X) X(:,3) > 0.109));   % boundary dofs, indices into q
[Kc, Mc, Vc] = fembd.craigBampton(K, M, ib, 20);       % u = T*Vc*[u_b; eta]
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
- **In the MATLAB Help browser:** once the toolbox folder is on the path, the documentation appears under *Supplemental Software*. Open a function's reference page with `doc fembd.Assembly`, `doc fembd.Body`, etc. The pages are prebuilt in `doc/html`; after editing `doc/src`, the tutorials or the examples, rebuild them with `cd doc; buildDocs`, which runs the tutorials and needs Gmsh.
- **Tutorials:** [`tutorials/`](tutorials/README.md), nine runnable scripts from a first model to joints, dynamics and Craig–Bampton reduction.
- **Manual:** [`docs/fembd.md`](docs/fembd.md) covers the workflow, formulation (elements, ties, elastic joints, rigid bodies, constraint elimination, Craig–Bampton), the API, recipes (loads, gravity, FRF, state space) and limitations.
- **In MATLAB:** `help fembd`, `help fembd.Assembly/tie`, etc.
- **Examples:** `examples/demoAssembly.m` (elastic parts plus a rigid payload: modes, Craig–Bampton, state space) `examples/demoParallelFlexure.m` (parallel leaf-spring guide between rigid bodies: 6×6 shuttle stiffness and modes, checked against beam theory) and `examples/demoSerialArm.m` (two-link arm: couplings and revolute joints with springs, reduced to joint coordinates).
- **Verification:** `tests/runTests.m`.

Requires MATLAB R2022a or later and Gmsh 4.
