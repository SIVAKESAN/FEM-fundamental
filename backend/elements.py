"""Compatible high-order 1D axial, 2D/3D truss and Euler–Bernoulli frame elements.

Internal hierarchical DOFs are fully retained in the global system; h-refined
members share ordinary physical nodes. No p-refinement is faked by submeshing.
"""
from __future__ import annotations
from functools import lru_cache
import numpy as np
import sympy as sp

T=sp.Symbol('t',real=True)

@lru_cache(None)
def axial_basis(p: int):
    # Endpoint-sharing hierarchical polynomial C0 basis: two linear endpoints
    # followed by p-1 internal bubbles. All interior modes vanish at endpoints.
    v=[1-T,T]
    for k in range(p-1):
        v.append(sp.expand(T*(1-T)*sp.legendre(k,2*T-1)))
    return tuple(map(sp.expand,v))

@lru_cache(None)
def bending_basis(p: int, length: float):
    # Cubic Hermite endpoint displacement and derivative interpolation.
    L=sp.Float(length,14)
    v=[1-3*T**2+2*T**3,
       L*(T-2*T**2+T**3),
       3*T**2-2*T**3,
       L*(-T**2+T**3)]
    for k in range(p-3):
        v.append(sp.expand(T**2*(1-T)**2*sp.legendre(k,2*T-1)))
    return tuple(map(sp.expand,v))

@lru_cache(None)
def eval_basis(p:int,length:float,bending:bool,derivative:int):
    shapes=bending_basis(p,length) if bending else axial_basis(p)
    return tuple(sp.lambdify(T,sp.diff(v,T,derivative),'numpy') for v in shapes)

def gauss_stiffness(basis, deriv, stiffness, L, p):
    n=len(deriv)
    # gauss n=p+2 integrates exact polynomials at p<=6; high precision.
    xs,ws=np.polynomial.legendre.leggauss(p+2)
    result=np.zeros((n,n))
    for z,w in zip(xs,ws):
        vals=np.array([float(f(float((z+1)/2))) for f in deriv],dtype=float)
        result+=float(w)/2*np.outer(vals,vals)
    return result * stiffness / (L ** (3 if basis=='bending' else 1))


def _dvals(p,L,bending,deriv,t):
    funcs=eval_basis(p,L,bending,deriv)
    return np.array([float(f(t)) for f in funcs])


def _section_local_stiffness(p,L,EA,EI=None,EIy=None,EIz=None,GJ=None,kind='bar1d'):
    """Build local element stiffness/consistent distributed member load in one basis.

    Returns k, f, basis description, endpoint-local map, and internals.
    local DOF ordering: axial endpoints u_i,u_j (and bubbles), other C1 groups.
    3D frame: endpoint vector [u,v,w,rx,ry,rz] at each end, then bubbles.
    """
    if kind in ('bar1d','truss2d','truss3d'):
        n=2+(p-1)
        K=np.zeros((n,n)); f=np.zeros(n)
        idx=list(range(n))
        K+=gauss_stiffness('axial',eval_basis(p,L,False,1),EA,L,p)
        return K,f,idx,n-2
    if kind=='beam1d':
        n=4+(p-3)
        K=np.zeros((n,n)); f=np.zeros(n)
        ii=list(range(n))
        K+=gauss_stiffness('bending',eval_basis(p,L,True,2),EI,L,p)
        return K,f,ii,n-4
    if kind=='frame2d':
        n=6+(p-1)+(p-3)
        K=np.zeros((n,n)); f=np.zeros(n)
        # local endpoints [u1, v1, rz1, u2, v2, rz2]
        axidx=[0,3]+list(range(6,6+p-1))
        v_idx=[1,2,4,5]+list(range(6+p-1,n))
        for family,idx,rigid,second in [('axial',axidx,EA,False),('bending',v_idx,EI,True)]:
            ks=gauss_stiffness(family,eval_basis(p,L,second,2 if second else 1),rigid,L,p)
            K[np.ix_(idx,idx)]+=ks
        return K,f,{'axial':axidx,'bending-y':v_idx},n-6
    if kind=='frame3d':
        n=12 + (p-1)*2 + (p-3)*2
        K=np.zeros((n,n)); f=np.zeros(n)
        # physical local [u,v,w,rx,ry,rz]*2
        # w Hermite slope = -ry -> factor - for local rotation y.
        axidx=[0,6]+list(range(12,12+p-1))
        txidx=[3,9]+list(range(12+p-1,12+2*(p-1)))
        vi=[1,5,7,11]+list(range(12+2*(p-1),12+2*(p-1)+(p-3)))
        wi=[2,4,8,10]+list(range(12+2*(p-1)+(p-3),n))
        for family,idx,rigid,second,signs in (
            ('axial',axidx,EA,False,None),('axial',txidx,GJ,False,None),
            ('bending',vi,EIz,True,None),
            ('bending',wi,EIy,True,np.array([1,-1,1,-1]+[1]*(p-3)))):
            ks=gauss_stiffness(family,eval_basis(p,L,second,2 if second else 1),rigid,L,p)
            if signs is not None:ks=ks*np.outer(signs,signs)
            K[np.ix_(idx,idx)]+=ks
        return K,f,{'axial':axidx,'torsion':txidx,'bending-local-y':vi,'bending-local-z':wi},n-12
    raise ValueError(f'Unsupported element {kind}')


def rotation_matrix_3d(a,b):
    ex=(np.array(b)-np.array(a));ex=ex/np.linalg.norm(ex)
    trial=np.array([0.0,0.0,1.0]) if abs(ex[2])<0.9 else np.array([0.0,1.0,0.0])
    ey=np.cross(trial,ex); ey=ey/np.linalg.norm(ey)
    ez=np.cross(ex,ey)
    return np.vstack((ex,ey,ez))


def make_element(kind,p,a,b,EA,EI,EIy,EIz,GJ,qx=0,qy=0,qz=0):
    a=np.array(a,dtype=float);b=np.array(b,dtype=float)
    direction=b-a
    L=float(np.linalg.norm(direction))
    k,f,idx,nprivate=_section_local_stiffness(p,round(L,12),EA,EI,EIy,EIz,GJ,kind)
    n=len(f)
    if kind in ('bar1d','truss2d','truss3d'):
        d={'bar1d':1,'truss2d':2,'truss3d':3}[kind]
        R=np.zeros((n,2*d+(p-1)))
        ex=direction[:d]/L
        R[0,0:d]=ex;R[1,d:2*d]=ex
        for i in range(p-1):R[i+2,2*d+i]=1
        Tglob=R
    elif kind=='beam1d':
        # Reversing member connectivity reverses the local x direction.
        # The local Hermite slope then changes sign relative to global rz.
        Tglob=np.eye(n)
        rotation_sign=1.0 if direction[0]>0 else -1.0
        Tglob[1,1]=rotation_sign
        Tglob[3,3]=rotation_sign
    elif kind=='frame2d':
        c,s=direction[0]/L,direction[1]/L
        Q=np.array([[c,s],[-s,c]])
        Tglob=np.eye(n)
        Tglob[np.ix_([0,1],[0,1])]=Q
        Tglob[np.ix_([3,4],[3,4])]=Q
    else:
        Q=rotation_matrix_3d(a,b)
        Tglob=np.eye(n)
        for at in (0,6):
            Tglob[at:at+3,at:at+3]=Q
            Tglob[at+3:at+6,at+3:at+6]=Q
    # Consistent distributed load vector by displacement interpolation.
    if kind in ('beam1d','frame2d','frame3d') and (qx or qy or qz):
        xs,ws=np.polynomial.legendre.leggauss(p+2)
        if kind=='beam1d':
            groups=[(list(range(4))+list(range(4,n)),qy,True)]
        elif kind=='frame2d':
            groups=[(idx['axial'],qx,False),(idx['bending-y'],qy,True)]
        else:
            groups=[(idx['axial'],qx,False),(idx['bending-local-y'],qy,True),(idx['bending-local-z'],qz,True)]
        for ii,q,bend in groups:
            if not q:continue
            signs=np.ones(len(ii))
            if kind=='frame3d' and bend and ii==idx['bending-local-z']:
                signs[1]=-1;signs[3]=-1
            for z,w in zip(xs,ws):
                vals=_dvals(p,round(L,12),bend,0,float((z+1)/2))*signs
                f[ii]+=q*L*float(w)/2*vals
    kg=Tglob.T@k@Tglob
    fg=Tglob.T@f
    return {'length':L,'direction':(direction/L).tolist(),'local':k,'global':kg,'T':Tglob,
            'f_local':f,'f_global':fg,'ndof_private':nprivate,'idx':idx,
            'axis':rotation_matrix_3d(a,b).tolist() if kind=='frame3d' else None}
