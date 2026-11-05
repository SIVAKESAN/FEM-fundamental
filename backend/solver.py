"""Small linear static FEM assembler supporting 1D/2D/3D beams, frames and trusses."""
from __future__ import annotations
import numpy as np
from scipy.linalg import solve as linsolve
from .models import Model, DOFS, FORCE_DOFS
from .elements import make_element

class FemError(ValueError):pass

def expand_model(model:Model):
    nodes=[n.model_dump() for n in model.nodes]
    elements=[]
    node_lookup={n['id']:n for n in nodes}
    next_id=max(node_lookup)+1
    for e in model.elements:
        a=node_lookup[e.i];b=node_lookup[e.j]
        ids=[e.i]
        for j in range(1,e.h):
            ratio=j/e.h
            new={'id':next_id,'x':a['x']+(b['x']-a['x'])*ratio,
                'y':a['y']+(b['y']-a['y'])*ratio,'z':a['z']+(b['z']-a['z'])*ratio,
                'fixed':[],'loads':{},
                'axis':[(b[k]-a[k])/max(sum((b[v]-a[v])**2 for v in ('x','y','z'))**0.5,1e-12) for k in ('x','y','z')]}
            nodes.append(new);node_lookup[next_id]=new
            ids.append(next_id);next_id+=1
        ids.append(e.j)
        for j in range(e.h):
            sub=e.model_dump()
            sub.update({'id':f'{e.id}.{j+1}' if e.h>1 else str(e.id),
                        'parent':e.id,'i':ids[j],'j':ids[j+1], 'h':1})
            elements.append(sub)
    return nodes,elements

def solve_model(model:Model,return_internal=False):
    kind=model.kind
    dofs=DOFS[kind];nn=len(dofs)
    nodes,els=expand_model(model)
    nmap={n['id']:n for n in nodes}
    original_ids={n.id for n in model.nodes}
    truss=kind in ('truss2d','truss3d')
    gmap={}
    description=[]
    for node in nodes:
        # A generated, collinear h-refinement truss joint has one axial displacement DOF.
        # Giving it unconstrained transverse DOFs would introduce fictitious mechanisms.
        local_keys=['ua'] if truss and node['id'] not in original_ids else dofs
        for key in local_keys:
            gmap[(node['id'],key)]=len(description)
            description.append(f"N{node['id']}.{key}")
    nphysical=len(description)
    g=nphysical
    members=[]
    for elem in els:
        a=nmap[elem['i']];b=nmap[elem['j']]
        data=make_element(kind,elem['p'],[a['x'],a['y'],a['z']],[b['x'],b['y'],b['z']],
                          elem['EA'],elem['EI'],elem['EIy'],elem['EIz'],elem['GJ'],
                          elem['qx'],elem['qy'],elem['qz'])
        endpoint_dofs=[]
        # Project a generated truss node's scalar axial movement back to global XYZ.
        # The projection also preserves the exact element-wise axial strain energy.
        if truss:
            dspace=len(dofs)
            columns=[]
            for n_id in (elem['i'],elem['j']):
                if n_id in original_ids:
                    endpoint_dofs.extend(gmap[(n_id,key)] for key in dofs)
                    columns.extend(np.eye(dspace))
                else:
                    endpoint_dofs.append(gmap[(n_id,'ua')])
                    columns.append(np.array(nmap[n_id]['axis'][:dspace]))
            # P maps compact global endpoints to all Cartesian coordinates used
            # by the ordinary two-ended truss transformation.
            count=len(endpoint_dofs);npriv=data['ndof_private']
            P=np.zeros((2*dspace+npriv,count+npriv))
            off=0;col=0
            for n_id in (elem['i'],elem['j']):
                if n_id in original_ids:
                    P[off:off+dspace,col:col+dspace]=np.eye(dspace)
                    col+=dspace
                else:
                    P[off:off+dspace,col]=nmap[n_id]['axis'][:dspace]
                    col+=1
                off+=dspace
            P[2*dspace:,count:]=np.eye(npriv)
            data['global']=P.T@data['global']@P
            data['f_global']=P.T@data['f_global']
            data['T']=data['T']@P
        else:
            endpoint_dofs=[gmap[(elem['i'],key)] for key in dofs]+[gmap[(elem['j'],key)] for key in dofs]
        edofs=endpoint_dofs
        for k in range(data['ndof_private']):
            edofs.append(g)
            description.append(f"E{elem['id']}.b{k+1}")
            g+=1
        if len(edofs)!=data['global'].shape[0]:
            raise FemError(f"Internal element DOF mismatch for element {elem['id']}")
        members.append({'element':elem,'data':data,'dofs':edofs})
    if g>550:raise FemError('The selected h/p refinement has too many DOFs (>550) for this educational solver. Reduce h or p.')
    K=np.zeros((g,g));F=np.zeros(g)
    for m in members:
        idx=m['dofs']
        K[np.ix_(idx,idx)]+=m['data']['global']
        F[idx]+=m['data']['f_global']
    fixed=[]
    for n in nodes:
        if truss and n['id'] not in original_ids:continue
        for d in n['fixed']:
            fixed.append(gmap[(n['id'],d)])
        for d in dofs:
            fkey=FORCE_DOFS[d]
            F[gmap[(n['id'],d)]]+=n['loads'].get(fkey,0.0)
    fixed=sorted(set(fixed))
    free=sorted(set(range(g))-set(fixed))
    if not fixed:
        raise FemError('No support restraints defined. Add fixed DOFs to eliminate rigid-body motion.')

    # H is the explicit boundary-condition matrix. Each row selects one
    # restrained global DOF: H @ u = prescribed_values (currently all zero).
    # Do not eliminate constrained DOFs: solve the Lagrange multiplier system.
    H=np.zeros((len(fixed),g),dtype=float)
    for row,global_dof in enumerate(fixed):
        H[row,global_dof]=1.0
    prescribed_values=np.zeros(len(fixed),dtype=float)

    # Check constraint independence, and catch disconnected free DOFs early.
    if np.linalg.matrix_rank(H)!=len(fixed):
        raise FemError('The boundary-condition matrix H contains dependent constraints.')
    diagonal=np.abs(np.diag(K))
    if any(diagonal[i]<1e-25 for i in free):
        raise FemError('A free DOF has no stiffness (mechanism). Add a support or connecting element.')

    # Numerical diagonal equilibration: let v=D*u, with D_ii=sqrt(|K_ii|).
    # For H containing unit rows, H*v = D_fixed*g, and the scaled multipliers
    # are mu=lambda/D_fixed. This is exactly the original saddle-point system
    # after a reversible change of variables, not a partitioned K_ff solve.
    positive=diagonal[diagonal>1e-25]
    reference=float(np.median(positive)) if positive.size else 1.0
    dscale=np.sqrt(np.where(diagonal>1e-25,diagonal,reference))
    Kscaled=K/dscale[:,None]/dscale[None,:]
    A_scaled=np.block([[Kscaled,H.T],[H,np.zeros((len(fixed),len(fixed)))]])
    rhs_scaled=np.concatenate((F/dscale,dscale[fixed]*prescribed_values))
    try:
        import warnings
        from scipy.linalg import LinAlgWarning
        with warnings.catch_warnings():
            warnings.simplefilter('error',LinAlgWarning)
            solution=linsolve(A_scaled,rhs_scaled,assume_a='sym',check_finite=True)
        if not np.all(np.isfinite(solution)):
            raise FemError('The augmented Lagrange multiplier matrix is singular or ill-conditioned.')
        relative_residual=np.linalg.norm(A_scaled@solution-rhs_scaled)/max(np.linalg.norm(rhs_scaled),1.0)
        if relative_residual>1e-7:
            raise FemError('The augmented Lagrange multiplier system has an excessive residual.')
    except Exception as exc:
        if isinstance(exc,FemError):raise
        raise FemError('Unable to solve the Lagrange multiplier system: the structure may be unstable or the constraints insufficient.') from exc

    u=solution[:g]/dscale
    multipliers=solution[g:]*dscale[fixed]
    R=K@u-F   # R = -H.T @ multipliers, by stationarity
    if not np.all(np.isfinite(u)):
        raise FemError('Numerical overflow in displacement solution')
    constraint_error=float(np.linalg.norm(H@u-prescribed_values))
    force_balance_err=float(np.linalg.norm((K@u-F)[free]))
    if constraint_error>1e-7:
        raise FemError('The Lagrange multiplier solution violates the displacement constraints.')
    node_results=[]
    for n in nodes:
        if truss and n['id'] not in original_ids:
            amplitude=float(u[gmap[(n['id'],'ua')]])
            dspl={d:amplitude*float(n['axis'][j]) for j,d in enumerate(dofs)}
            reactions={}
        else:
            dspl={d:float(u[gmap[(n['id'],d)]]) for d in dofs}
            reactions={FORCE_DOFS[d]:float(R[gmap[(n['id'],d)]]) for d in dofs if d in n['fixed']}
        node_results.append({'id':n['id'],'x':n['x'],'y':n['y'],'z':n['z'],
                             'displacements':dspl,'reactions':reactions,'fixed':n['fixed'],
                             'generated':n['id'] not in original_ids})
    member_results=[]
    for m in members:
        e=m['element'];d=m['data'];ed=u[m['dofs']]
        local=d['T']@ed
        forces=d['local']@local-d['f_local']
        member_results.append({'id':e['id'],'parent':e['parent'],'i':e['i'],'j':e['j'],
                               'p':e['p'],'L':d['length'],'local_dofs':_finite(local),
                               'end_forces_local':_finite(forces[:2] if truss or kind=='bar1d' else forces[:nn*2])})
    result={
        'kind':kind,'dofs':description,'n_physical_dofs':nphysical,'n_total_dofs':g,
        'n_restrained_dofs':len(fixed),'free_indices':free,'fixed_indices':fixed,
        'nodes':node_results,'elements':member_results,
        'F':_finite(F),'K':_matrix(K),
        'H':_matrix(H),'g':_finite(prescribed_values),'lambda':_finite(multipliers),
        'solver_method':'lagrange_multiplier','u':_finite(u),'R':_finite(R),
        'summary':{'max_translation_m':max((abs(float(n['displacements'][d])) for n in node_results for d in dofs if d in ('ux','uy','uz')),default=0),
                   'max_rotation_rad':max((abs(float(n['displacements'][d])) for n in node_results for d in dofs if d.startswith('r')),default=0),
                   'free_equilibrium_residual_N':force_balance_err,
                   'constraint_residual':constraint_error,
                   'lagrange_equilibrium_residual':float(np.linalg.norm(K@u+H.T@multipliers-F)),
                   'elements_after_h_refinement':len(els),
                   'status':'Solved — linear elastic, static, small deformation'}
    }
    if return_internal:result['_internal']={'members':members,'nodes':nodes,'K':K,'F':F,'u':u,'R':R,'free':free,'fixed':fixed,
                                           'H':H,'g':prescribed_values,'lambda':multipliers}
    return result

def _finite(v):return [float(x) for x in v]
def _matrix(M):return [[float(x) for x in row] for row in M]
