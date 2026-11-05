"""Readable, model-specific FEM solution manual from the *actual* LM solver.

Order: interpolation -> Ke -> T -> KG -> H -> stationary Lagrangian ->
augmented matrix -> solve -> reactions -> member forces.
All numerical matrices are extracted from the same assembly used by the solver.
"""
from __future__ import annotations
import numpy as np
import sympy as sp
from sympy.printing.mathml import mathml
from .models import Model, DOFS
from .elements import T, axial_basis, bending_basis
from .solver import solve_model, FemError


def _mml(expr):
    return '<math xmlns="http://www.w3.org/1998/Math/MathML" display="block">'+mathml(expr,printer='presentation')+'</math>'


def _eq(name, value):
    return _mml(sp.Eq(sp.Symbol(name),value,evaluate=False))


def _formula(inner):
    """Native presentation MathML for canonical matrix formulas (not numeric data)."""
    return '<math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mrow>'+inner+'</mrow></math>'


def _var(name, sub=None, sup=None):
    base=f'<mi>{name}</mi>'
    if sub is not None: base=f'<msub>{base}<mi>{sub}</mi></msub>'
    if sup is not None: base=f'<msup>{base}<mi>{sup}</mi></msup>'
    return base


def _standard_equations():
    d=_var('d'); dt=_var('d',sup='T'); kg=_var('K',sub='G'); f=_var('F')
    h=_var('H'); ht=_var('H',sup='T'); lam=_var('λ'); lamt=_var('λ',sup='T'); g=_var('g')
    pi=_formula(_var('Π')+'<mo>=</mo><mfrac><mn>1</mn><mn>2</mn></mfrac>'+dt+kg+d+
                '<mo>−</mo>'+dt+f+'<mo>+</mo>'+lamt+'<mo>(</mo>'+h+d+'<mo>−</mo>'+g+'<mo>)</mo>')
    equilibrium=_formula(kg+d+'<mo>+</mo>'+ht+lam+'<mo>=</mo>'+f)
    constraints=_formula(h+d+'<mo>=</mo>'+g)
    reactions=_formula(_var('R')+'<mo>=</mo><mo>−</mo>'+ht+lam)
    return pi,equilibrium,constraints,reactions


def matrix_math(name, M, max_dim=26):
    """Presentation MathML of a real numerical matrix, or None for large systems."""
    mat=np.asarray(M,dtype=float)
    if mat.ndim==1: mat=mat.reshape((-1,1))
    if max(mat.shape)>max_dim: return None
    cutoff=max(1e-12,1e-14*float(np.max(np.abs(mat)))) if mat.size else 1e-12
    def tidy(x):
        x=float(x)
        if abs(x)<cutoff: return sp.Integer(0)
        if abs(x-round(x))<1e-12: return sp.Integer(round(x))
        if abs(x)>=1e6 or 0<abs(x)<1e-4: return sp.Float(x,5)
        return sp.Float(x,7)
    mm=sp.Matrix(mat.shape[0],mat.shape[1],lambda i,j:tidy(mat[i,j]))
    return _eq(name,mm)


def derive(model:Model, selected_id:str):
    res=solve_model(model,return_internal=True)
    intr=res.pop('_internal')
    member=next((m for m in intr['members'] if str(m['element']['id'])==selected_id),None)
    if member is None:
        member=next((m for m in intr['members'] if str(m['element']['parent'])==selected_id),None)
    if member is None: raise FemError('Select a valid element ID')

    e=member['element'];kind=model.kind;data=member['data'];p=e['p'];L=data['length']
    bending=kind in ('beam1d','frame2d','frame3d')
    n=len(res['dofs']); nc=len(intr['fixed'])
    steps=[]
    def step(title,body,equations=(),note=''):
        steps.append({'title':title,'body':body,
                      'equations':[x for x in equations if x], 'note':note})

    # 1. The actual compatible basis used by the numerical element.
    sf=[]
    if kind!='beam1d':
        for i,v in enumerate(axial_basis(p)):
            sf.append(_eq(f'N_{{{i+1}}}(t)',v))
    hf=[]
    if bending:
        for i,v in enumerate(bending_basis(p,round(L,12))):
            hf.append(_eq(f'h_{{{i+1}}}(t)',v))
    note=(f"Element {e['id']} ({e['i']} → {e['j']}), length {L:.6g} m, p={p}. "
          f"Physical DOFs per node: {', '.join(DOFS[kind])}. "
          "The reference coordinate is t ∈ [0,1], with x_local=L·t and Jacobian J=L. "
          "Hierarchical interior p-modes are real unknowns in the global system.")
    step('1. Interpolation and shape functions',
         'Start from the displacement approximation. Axial and torsional modes use C⁰ interpolation; Euler–Bernoulli bending uses C¹ Hermite interpolation. Each equation below is the actual basis used by this element.',
         [_eq('J',sp.Float(L,8))]+sf+hf,note)

    # 2. Derive Ke by differentiating the actual basis, then show its value.
    Lsym=sp.Symbol('L',positive=True)
    EAsym=sp.Symbol('EA',positive=True)
    EIsym=sp.Symbol('EI',positive=True)
    EIysym=sp.Symbol('EI_y',positive=True)
    EIzsym=sp.Symbol('EI_z',positive=True)
    GJsym=sp.Symbol('GJ',positive=True)
    Ni=sp.Function('N_i')(T); Nj=sp.Function('N_j')(T)
    Hi=sp.Function('H_i')(T); Hj=sp.Function('H_j')(T)
    eqs=[]
    if kind!='beam1d':
        eqs.append(_eq('B_{axial}(t)',sp.Matrix([sp.diff(v,T)/sp.Float(L,12) for v in axial_basis(p)]).T))
        eqs.append(_mml(sp.Eq(sp.Symbol('k_{ij}^{axial}'),
            EAsym/Lsym*sp.Integral(sp.diff(Ni,T)*sp.diff(Nj,T),(T,0,1)),evaluate=False)))
        if p==1:
            eqs.append(_eq('K_{e,axial}', EAsym/Lsym*sp.Matrix([[1,-1],[-1,1]])))
    if bending:
        basis=bending_basis(p,round(L,12))
        eqs.append(_eq('B_{bending}(t)',sp.Matrix([sp.diff(v,T,2)/sp.Float(L,12)**2 for v in basis]).T))
        eqs.append(_mml(sp.Eq(sp.Symbol('k_{ij}^{bending}'),
            EIsym/Lsym**3*sp.Integral(sp.diff(Hi,T,2)*sp.diff(Hj,T,2),(T,0,1)),evaluate=False)))
        if kind=='frame3d':
            eqs.append(_mml(sp.Eq(sp.Symbol('k_{ij}^{bending,y}'),
                EIysym/Lsym**3*sp.Integral(sp.diff(Hi,T,2)*sp.diff(Hj,T,2),(T,0,1)),evaluate=False)))
            eqs.append(_mml(sp.Eq(sp.Symbol('k_{ij}^{bending,z}'),
                EIzsym/Lsym**3*sp.Integral(sp.diff(Hi,T,2)*sp.diff(Hj,T,2),(T,0,1)),evaluate=False)))
            eqs.append(_mml(sp.Eq(sp.Symbol('k_{ij}^{torsion}'),
                GJsym/Lsym*sp.Integral(sp.diff(Ni,T)*sp.diff(Nj,T),(T,0,1)),evaluate=False)))
        if p==3:
            cl=sp.Matrix([[12,6*Lsym,-12,6*Lsym],
                          [6*Lsym,4*Lsym**2,-6*Lsym,2*Lsym**2],
                          [-12,-6*Lsym,12,-6*Lsym],
                          [6*Lsym,2*Lsym**2,-6*Lsym,4*Lsym**2]])
            eqs.append(_eq('K_{e,bending}',EIsym/Lsym**3*cl))
    eqs.append(matrix_math('K_{e,local}',data['local']))
    eqs.append(matrix_math('f_{e,local}',data['f_local']))
    step('2. Derive element stiffness Kₑ',
         'Differentiate the shape functions to form the strain–displacement operator B. Integrate strain energy to obtain each element stiffness entry. The final numerical matrix is the actual Gaussian-integrated local Kₑ. The local consistent load vector is also shown.',
         eqs,
         'For 3D frames, the actual local Kₑ additionally combines EA, GJ, EIy and EIz in separate axial, torsional and bending blocks. Bending about local y uses the negative local ry slope convention. Stiffness values are in SI mixed force/displacement and moment/rotation units.')

    # 3. The actual member coordinate transformation.
    eqs=[_formula(_var('K',sub='e,G')+'<mo>=</mo>'+_var('T',sup='T')+_var('K',sub='e,local')+_var('T')),
         matrix_math('T',data['T']),matrix_math('K_{e,G}',data['global']),
         matrix_math('f_{e,G}',data['f_global'])]
    step('3. Local-to-global transformation',
         'For trusses and inclined frames, resolve local member displacement components into global nodal components. Transform both stiffness and consistent loads. For a 1D member, T is the identity or axial-direction projection.',
         eqs,
         'u_local = T·u_element; K_e,G = Tᵀ·K_e,local·T; f_e,G = Tᵀ·f_e,local. The displayed T is the exact transformation used by the solver, including hierarchical internal DOFs.')

    # 4. Assembly with explicit DOF indices and real matrices.
    element_indices=[j+1 for j in member['dofs']]
    mapping=', '.join(f'{j+1}={res["dofs"][j]}' for j in member['dofs'])
    step('4. Assemble global stiffness Kᴳ and force F',
         'Use each element connectivity mapping to add its transformed stiffness entries into the global Kᴳ. Contributions at shared DOFs are summed. Assemble nodal forces and equivalent element loads into F.',
         [_eq('I_e',sp.Matrix(element_indices)),
          matrix_math('K_G',intr['K']),matrix_math('F',intr['F'])],
         f"Selected element assembly indices (1-based): {mapping}. "
         f"The global stiffness matrix is {n}×{n}. Full numerical matrices are displayed up to 26 rows/columns; larger ones are still assembled and solved.")

    # 5. Constraint matrix H, not a reduced free-DOF matrix.
    fixed_labels=[res['dofs'][i] for i in intr['fixed']]
    step('5. Construct boundary-condition matrix H',
         'For every restrained displacement DOF, place a 1 in the corresponding column of H and zeros elsewhere. The support conditions are H·d=g. The current model supports zero prescribed displacements, so g=0.',
         [_mml(sp.Eq(sp.Symbol('H')*sp.Symbol('d'),sp.Symbol('g'),evaluate=False)),
          matrix_math('H',intr['H']),matrix_math('g',intr['g'])],
         f"H has {nc} rows (one per independent constraint) and {n} columns (one per global DOF). "
         f"Constraint row order: {', '.join(f'{k+1}: {label}=0' for k,label in enumerate(fixed_labels))}. "
         'The Hermite bending shape functions are not this boundary-condition matrix H.')

    # 6. Lagrange multiplier stationary potential.
    pi_formula,eq_formula,constraint_formula,reaction_formula=_standard_equations()
    step('6. Form the Lagrangian and stationary equations',
         'Enforce supports without deleting the restrained DOFs. Add one multiplier for each independent row of H to the total potential energy. Differentiate with respect to d and λ.',
         [pi_formula,eq_formula,constraint_formula],
         'With the sign convention Π=½dᵀKᴳd−dᵀF+λᵀ(Hd−g), the physical support reaction vector is R=−Hᵀλ.')

    # 7. Explicit augmented matrix in *unscaled physical units*.
    K=intr['K']; H=intr['H']; F=intr['F']; gv=intr['g']
    A=np.block([[K,H.T],[H,np.zeros((nc,nc))]])
    rhs=np.concatenate((F,gv))
    symbol_A=sp.Matrix([[sp.Symbol('K_G'),sp.Symbol('H^T')],
                        [sp.Symbol('H'),sp.Integer(0)]])
    symbol_x=sp.Matrix([sp.Symbol('d'),sp.Symbol('λ')])
    symbol_rhs=sp.Matrix([sp.Symbol('F'),sp.Symbol('g')])
    step('7. Assemble the Lagrange multiplier system',
         'Combine the original global stiffness matrix and the H matrix into one augmented saddle-point system. This is the matrix actually solved by the Python backend, with numerical diagonal equilibration for stability.',
         [_mml(sp.Eq(sp.MatMul(symbol_A,symbol_x,evaluate=False),symbol_rhs,evaluate=False)),
          matrix_math('A_{LM}',A),matrix_math('b_{LM}',rhs)],
         f"Augmented matrix size: ({n}+{nc})×({n}+{nc}) = {n+nc}×{n+nc}. "
         'The displayed A_LM and b_LM are in original physical units. Internal diagonal scaling does not change the mathematical system or solution.')

    # 8. Solution: displacements AND multipliers.
    step('8. Solve for nodal displacements and multipliers',
         'Solve the augmented equations simultaneously. The first entries are all displacement DOFs (including internal p-refinement modes); the final entries are Lagrange multipliers corresponding to H rows.',
         [matrix_math('d',intr['u']),matrix_math('λ',intr['lambda']),
          matrix_math('H d - g',H@intr['u']-gv)],
         f"Constraint residual: {res['summary']['constraint_residual']:.3g}; "
         f"augmented equilibrium residual: {res['summary']['lagrange_equilibrium_residual']:.3g}. "
         'For fixed-zero constraints, the selected displacement components should be zero within numerical precision.')

    # 9. Reactions with explicit sign convention.
    reactions_from_lambda=-H.T@intr['lambda']
    step('9. Recover support reactions and verify equilibrium',
         'Calculate the reaction vector from the original full stiffness matrix and force vector. With the chosen Lagrangian sign, it must equal −Hᵀλ. Only constrained DOFs have support reactions.',
         [_formula(_var('R')+'<mo>=</mo>'+_var('K',sub='G')+_var('d')+'<mo>−</mo>'+_var('F')),
          reaction_formula,
          matrix_math('R',intr['R']),
          matrix_math('R_{LM}',reactions_from_lambda)],
         f"Maximum absolute reaction-vector difference: {np.max(np.abs(intr['R']-reactions_from_lambda)):.4g}. "
         f"Free-DOF equilibrium residual norm: {res['summary']['free_equilibrium_residual_N']:.4g}.")

    # 10. Selected element internal forces.
    d_e=intr['u'][member['dofs']]
    d_local=data['T']@d_e
    q_local=data['local']@d_local-data['f_local']
    step('10. Recover selected element forces',
         'Extract this element’s solved global displacement vector, transform to local coordinates and recover its local end-force vector. For higher-order elements, internal modal forces are also present and should satisfy modal equilibrium.',
         [matrix_math('d_e',d_e),matrix_math('d_{e,local}',d_local),
          _mml(sp.Eq(sp.Symbol('q_{e,local}'),sp.Symbol('K_{e,local}d_{e,local}-f_{e,local}'),evaluate=False)),
          matrix_math('q_{e,local}',q_local)],
         'For bars and trusses, the first two local entries correspond to axial endpoint forces. For beam/frame elements, the endpoint entries are ordered by the physical local nodal DOFs.')

    return {'element':e['id'],'kind':kind,'steps':steps,
            'solver_method':'lagrange_multiplier',
            'element_dofs':[res['dofs'][j] for j in member['dofs']],
            'constraint_dofs':fixed_labels,
            'selected_geometry':{'length_m':L,'p':p,
              'h':len([m for m in intr['members'] if m['element']['parent']==e['parent']])}}
