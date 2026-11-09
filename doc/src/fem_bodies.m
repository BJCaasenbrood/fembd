%% Bodies, Meshing and Materials
% A body is one part of the assembly: a mesh and a material. It can be
% elastic or rigid.
%
%% Creating a body
%
%   b = fem.Body(src, mat, Name=Value)
%
% |src| is an STL file (meshed with Gmsh), an ASCII Gmsh .msh file (version
% 2.2 or 4.1), or a struct with fields |X| (n x 3 nodes) and |T| (ne x 4
% or ne x 10 tetrahedra, Gmsh node order).
%
% <html>
% <table border="1" cellpadding="4" style="border-collapse:collapse">
% <tr><th>Option</th><th>Default</th><th>Description</th></tr>
% <tr><td>Name</td><td>file name</td><td>Name used to refer to the body in an assembly</td></tr>
% <tr><td>Scale</td><td>1</td><td>Coordinate factor applied after meshing, e.g. 1e-3 for mm to m</td></tr>
% <tr><td>Size</td><td>Gmsh default</td><td>Maximum element size, in file units</td></tr>
% <tr><td>Order</td><td>2</td><td>2 = quadratic Tet10, 1 = linear Tet4</td></tr>
% <tr><td>Remesh</td><td>true</td><td>Remesh the STL surface (false: keep the STL triangles)</td></tr>
% <tr><td>FeatureAngle</td><td>40</td><td>Angle (deg) used to split the STL into smooth patches</td></tr>
% <tr><td>Rigid</td><td>false</td><td>Make the body rigid (6 dofs at RefPoint)</td></tr>
% <tr><td>RefPoint</td><td>centroid</td><td>Reference point of a rigid body</td></tr>
% <tr><td>Damping</td><td>[0 0]</td><td>Rayleigh damping [alpha beta] of an elastic body: D = alpha*M + beta*K</td></tr>
% </table>
% </html>
%
%% Element order
% Use the default Tet10 elements. Linear Tet4 elements are far too stiff in
% bending: on a cantilever, Tet4 is still about 9% off on a fine mesh where
% Tet10 is within 0.2% on a coarse one (see
% <tut02_meshing.html Tutorial 2>).
%
% Tet10 edge nodes are placed on the curved surface reconstructed from the
% STL. If that inverts elements (coarse meshes on small radii), meshSTL
% remeshes with straight-sided Tet10 and warns; a smaller |Size| keeps the
% curved geometry.
%
%% Materials
% A material is a struct with a 6 x 6 stiffness |D| in Voigt notation
% (order xx yy zz yz xz xy, engineering shear strains) and a density
% |rho|. |fem.material(E, nu, rho)| gives an isotropic material; any
% symmetric positive-definite |D| works, e.g. orthotropic composites (see
% <tut03_materials.html Tutorial 3>).
%
%   mat = fem.material(210e9, 0.3, 7850);
%   cfrp = struct('D', D, 'rho', 1600);
%
%% Mass matrix
% |b.matrices()| returns the body's K and consistent M. With
% |Lumped = true| (also in |a.system|) M is diagonal, using HRZ lumping:
% each element's mass diagonal is scaled to keep the element mass.
%
%% Rigid bodies
% With |Rigid = true| the body has no stiffness and 6 dofs
% [ux uy uz rx ry rz] at its reference point; all its nodes follow
% u = u_ref + r x (x - x_ref). The mesh only provides mass, inertia and
% the surface for ties and plots, so a coarse Tet4 mesh is enough.
%
% Rigid bodies do not need a volume mesh at all:
%
%   b = fem.Body.rigidSurface('lid.stl', 7850, Scale = 1e-3);   % STL surface only
%   b = fem.Body.pointMass(m, c, J, Name = 'motor');            % mass properties only
%
% |rigidSurface| computes mass, centre of mass and inertia exactly from
% the closed STL surface. |pointMass| takes mass |m|, centre of mass |c|
% and inertia tensor |J| (about |c|), e.g. from a CAD report; it has no
% surface and is connected with a joint. See
% <tut06_rigidBodies.html Tutorial 6>.
%
% To use measured or CAD mass properties instead, call |setMass| on any
% rigid body. The inertia can be given about another point, in principal
% axes, or in mm, and is checked for validity:
%
%   b.setMass(1.25);                                        % scale: same distribution
%   b.setMass(1.25, CenterOfMass = c, Inertia = J);         % J about c, global axes
%   b.setMass(1.25, CenterOfMass = c, Inertia = Jo, InertiaPoint = o);   % J about point o
%   b.setMass(1.25, CenterOfMass = c, Inertia = [J1 J2 J3], Axes = R);   % principal axes
%   b.setMass(1.25, CenterOfMass = cmm, Inertia = Jmm, LengthScale = 1e-3);   % mm, kg mm^2
%
%% Placing parts
% |b.transform(R, t)| moves a body rigidly, x -> R*x + t.
%
%% See also
% <matlab:doc('fem.Body') fem.Body>, <matlab:doc('fem.meshSTL') fem.meshSTL>,
% <matlab:doc('fem.readMsh') fem.readMsh>, <matlab:doc('fem.material') fem.material>
