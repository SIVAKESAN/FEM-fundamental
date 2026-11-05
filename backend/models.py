"""Typed model validation for the educational linear FEM application.

Base SI units: metres, newtons, radians, N, N*m, EA in N,
EI and GJ in N*m^2. All constraints are zero-valued displacement DOFs.
"""
from __future__ import annotations
from typing import Literal
from pydantic import BaseModel, Field, model_validator

Kind = Literal['bar1d', 'truss2d', 'truss3d', 'beam1d', 'frame2d', 'frame3d']
DOFS = {
    'bar1d': ['ux'],
    'truss2d': ['ux', 'uy'],
    'truss3d': ['ux', 'uy', 'uz'],
    'beam1d': ['uy', 'rz'],
    'frame2d': ['ux', 'uy', 'rz'],
    'frame3d': ['ux', 'uy', 'uz', 'rx', 'ry', 'rz'],
}
FORCE_DOFS = {'ux':'fx','uy':'fy','uz':'fz','rx':'mx','ry':'my','rz':'mz'}

class Node(BaseModel):
    id: int = Field(ge=1)
    x: float = 0
    y: float = 0
    z: float = 0
    fixed: list[str] = Field(default_factory=list)
    loads: dict[str, float] = Field(default_factory=dict)

class Element(BaseModel):
    id: int = Field(ge=1)
    i: int
    j: int
    EA: float = Field(default=2e7, gt=0)
    EI: float = Field(default=2e5, gt=0)
    EIy: float = Field(default=2e5, gt=0)
    EIz: float = Field(default=2e5, gt=0)
    GJ: float = Field(default=1e5, gt=0)
    p: int = Field(default=3, ge=1, le=6)
    h: int = Field(default=1, ge=1, le=8)
    qx: float = 0.0
    qy: float = 0.0
    qz: float = 0.0

class Model(BaseModel):
    kind: Kind = 'bar1d'
    nodes: list[Node]
    elements: list[Element]

    @model_validator(mode='after')
    def validate_model(self):
        if not self.nodes or not self.elements:
            raise ValueError('At least one node and one element are required')
        node_ids = [n.id for n in self.nodes]
        if len(node_ids) != len(set(node_ids)):
            raise ValueError('Node IDs must be unique')
        element_ids = [e.id for e in self.elements]
        if len(element_ids) != len(set(element_ids)):
            raise ValueError('Element IDs must be unique')
        allowed = set(DOFS[self.kind]); loads_allowed = {FORCE_DOFS[k] for k in allowed}
        for n in self.nodes:
            if set(n.fixed)-allowed:
                raise ValueError(f'Node {n.id}: invalid restrained DOF for {self.kind}: {set(n.fixed)-allowed}')
            if set(n.loads)-loads_allowed:
                raise ValueError(f'Node {n.id}: invalid load component for {self.kind}: {set(n.loads)-loads_allowed}')
        coords = {n.id: (n.x,n.y,n.z) for n in self.nodes}
        for e in self.elements:
            if e.i not in coords or e.j not in coords:
                raise ValueError(f'Element {e.id} references nonexistent node')
            if e.i == e.j:
                raise ValueError(f'Element {e.id} has identical start and end nodes')
            xa,xb=coords[e.i],coords[e.j]
            if self.kind in ('bar1d','beam1d'):
                d=(xb[0]-xa[0])**2
                if abs(xb[1]-xa[1]) > 1e-9 or abs(xb[2]-xa[2]) > 1e-9:
                    raise ValueError(f'Element {e.id}: 1D bars/beams must align with global x')
            elif self.kind in ('truss2d','frame2d'):
                d=sum((xb[j]-xa[j])**2 for j in (0,1))
                if abs(xb[2]-xa[2])>1e-9:
                    raise ValueError(f'Element {e.id}: 2D members must lie in XY plane')
            else:
                d=sum((xb[j]-xa[j])**2 for j in range(3))
            if d <= 1e-18:
                raise ValueError(f'Element {e.id} has zero or near-zero length')
            if self.kind in ('beam1d','frame2d','frame3d') and e.p < 3:
                raise ValueError(f'Element {e.id}: Euler–Bernoulli bending requires p ≥ 3 (cubic Hermite)')
            if self.kind in ('bar1d','truss2d','truss3d') and (e.qx or e.qy or e.qz):
                raise ValueError('For this version, distributed member loading is supported only for beams/frames')
        if len(self.nodes)>80 or sum(e.h for e in self.elements)>160:
            raise ValueError('Model size exceeds the educational demo limit (80 input nodes, 160 sub-elements)')
        return self
