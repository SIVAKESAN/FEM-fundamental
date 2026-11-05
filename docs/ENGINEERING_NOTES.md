# Engineering and mathematical notes

**Units:** global translations [m], rotations [rad], nodal forces [N], nodal moments [N·m], EA [N], EI and GJ [N·m²].

## Axial bar and truss

Axial strain energy: `U = (1/2) integral EA (du/dx)^2 dx`.

Let `t=x/L` and use `u = N_i(t) a_i`. The axial stiffness is

`K_local(i,j) = (EA/L) integral_0^1 N'_i(t) N'_j(t) dt`.

For p=1 the classic local stiffness is

`K_local = (EA/L) [[1,-1],[-1,1]]`.

A 2D/3D truss endpoint axial DOF is the component of its nodal displacement in the unit member direction, e.g. `u_local,1 = c_x ux_1 + c_y uy_1 + c_z uz_1`. Thus `K_global = T^T K_local T`. Internal p-refinement DOFs are axial, not artificial bending degrees of freedom.

## Euler–Bernoulli beam

The compatible interpolation is C1 across member endpoints. The cubic Hermite base has displacement and slope as nodal DOFs. Bending strain energy is

`U_b = (1/2) integral EI (d²w/dx²)² dx`, and

`K_b(i,j) = EI/L³ integral_0^1 H''_i(t) H''_j(t) dt`.

For p>3, internal polynomials of the form `t²(1-t)² P_k(2t-1)` provide C1-compatible higher-order bending modes. These vanish in both value and first derivative at element ends, preserving nodal continuity. The p-refinement is genuine, not an alias for h-refinement.

## Plane / space frame

- Plane frame nodal DOFs `[ux, uy, rz]`.
- Space frame nodal DOFs `[ux,uy,uz,rx,ry,rz]`.
- Axial stiffness `EA` and torsion `GJ` are first-derivative energies.
- Bending stiffness about local z is `EIz`, about local y is `EIy`.
- The 3D local basis is a right-handed orthonormal coordinate system, dynamically constructed from member endpoints. Angular DOFs rotate with the same 3D basis as translational DOFs.

## Refinement

- **h-refinement:** split a physical member into h collinear sub-members with shared physical displacements, rotations (for beams/frames) or a single axial joint DOF for intermediate truss points.
- **p-refinement:** increase compatible shape-function polynomial order, introducing independent hierarchical internal element unknowns. Both refinements are solved in the global stiffness system.

## Boundary conditions and the Lagrange multiplier method

**The solver does not use the partitioned free-DOF matrix `Kff` as its solution method.**

After deriving element matrices and assembling `K_G` and `F`, create a boundary-condition matrix `H`. For every restrained global displacement DOF `j`, insert a constraint row with `H[i,j]=1` and zeros elsewhere. For the current zero-support-settlement model, `g=0` and `H d = g`.

Introduce the Lagrangian:

`Pi(d,lambda) = (1/2) d^T K_G d - d^T F + lambda^T (H d - g)`.

Stationarity with respect to displacements and multipliers gives:

`K_G d + H^T lambda = F` and `H d = g`.

Hence the augmented saddle-point system is

```text
[ K_G   H^T ] [ d      ] = [ F ]
[ H      0  ] [ lambda ]   [ g ]
```

**Reaction sign convention:** the full global reaction residual is

`R = K_G d - F = -H^T lambda`.

For an identity-row support constraint, the multiplier `lambda_i` is the negative of the corresponding physical reaction. This is verified numerically for every solved model.

The implementation uses reversible diagonal scaling of the *augmented* system to improve numerical conditioning, not a reduced `Kff` solve. The displayed matrices are in original physical units. Constraints are checked for independence, mechanisms are rejected, and the constraint and equilibrium residuals are reported.

Currently the UI supports zero-prescribed restraints only. Nonzero support movements, general multipoint constraints, releases, nonlinear behaviour, mass matrices and geometric second-order effects are outside this version.

## Formula display

The **SHOW COMPLETE SOLUTION** panel renders SymPy-generated **presentation MathML** as real equations in current Chrome/Edge, including the numerical stiffness and transformation matrices. Matrix displays above 26 DOFs are omitted for readability but remain in the solver results.
