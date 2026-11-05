import copy
import numpy as np
import pytest
from fastapi.testclient import TestClient
from backend.models import Model
from backend.solver import solve_model, FemError
from backend.examples import EXAMPLES
from backend.app import app
from backend.derivations import derive

@pytest.mark.parametrize('kind',list(EXAMPLES))
def test_example_solves(kind):
    r=solve_model(Model.model_validate(EXAMPLES[kind]))
    assert r['summary']['free_equilibrium_residual_N']<0.00001
    assert np.all(np.isfinite(r['u']))
    assert len(r['nodes'])==len(EXAMPLES[kind]['nodes'])

@pytest.mark.parametrize('kind',list(EXAMPLES))
def test_h_refinement_is_stable(kind):
    base=copy.deepcopy(EXAMPLES[kind]);fine=copy.deepcopy(base)
    for e in fine['elements']:e['h']=2
    r0=solve_model(Model.model_validate(base))
    r1=solve_model(Model.model_validate(fine))
    for n0,n1 in zip(r0['nodes'],r1['nodes'][:len(base['nodes'])]):
        for d in n0['displacements']:
            assert n1['displacements'][d]==pytest.approx(n0['displacements'][d],abs=1e-7,rel=1e-6)

@pytest.mark.parametrize('kind',list(EXAMPLES))
def test_p_refinement_is_stable(kind):
    low=copy.deepcopy(EXAMPLES[kind]);high=copy.deepcopy(low)
    for e in low['elements']:e['p']=3 if kind in ('beam1d','frame2d','frame3d') else 1
    for e in high['elements']:e['p']=4
    r0=solve_model(Model.model_validate(low));r1=solve_model(Model.model_validate(high))
    for n0,n1 in zip(r0['nodes'],r1['nodes']):
        for d in n0['displacements']:
            assert n1['displacements'][d]==pytest.approx(n0['displacements'][d],abs=1e-7,rel=1e-6)

def test_indeterminate_axial_bar():
    # 4 m fixed-fixed bar, central concentrated load; equal end EA/L stiffness.
    m=copy.deepcopy(EXAMPLES['bar1d'])
    res=solve_model(Model.model_validate(m))
    center=res['nodes'][1]
    assert center['displacements']['ux']==pytest.approx(10000/(2*(2e7/2)))
    assert res['nodes'][0]['reactions']['fx']==pytest.approx(-5000)
    assert res['nodes'][2]['reactions']['fx']==pytest.approx(-5000)

def test_cantilever_beam_tip_load():
    EI=2e5;P=1000.;L=4.
    m={'kind':'beam1d','nodes':[{'id':1,'x':0,'fixed':['uy','rz']},
        {'id':2,'x':L,'loads':{'fy':-P}}],
        'elements':[{'id':1,'i':1,'j':2,'EI':EI,'p':3}]}
    r=solve_model(Model.model_validate(m))
    assert r['nodes'][1]['displacements']['uy']==pytest.approx(-P*L**3/(3*EI))
    assert r['nodes'][1]['displacements']['rz']==pytest.approx(-P*L**2/(2*EI))

def test_simply_supported_beam_center_load():
    L=6.;EI=2e5;P=1500
    m={'kind':'beam1d','nodes':[{'id':1,'x':0,'fixed':['uy']},
        {'id':2,'x':L/2,'loads':{'fy':-P}},
        {'id':3,'x':L,'fixed':['uy']}],
        'elements':[{'id':1,'i':1,'j':2,'EI':EI,'p':3},
                    {'id':2,'i':2,'j':3,'EI':EI,'p':3}]}
    r=solve_model(Model.model_validate(m))
    assert r['nodes'][1]['displacements']['uy']==pytest.approx(-P*L**3/(48*EI))
    assert r['nodes'][0]['reactions']['fy']==pytest.approx(P/2)
    assert r['nodes'][2]['reactions']['fy']==pytest.approx(P/2)

def test_cantilever_frame2d_tip_force():
    L=4.;EI=2e5;P=1000.
    m={'kind':'frame2d','nodes':[{'id':1,'x':0,'y':0,'fixed':['ux','uy','rz']},
        {'id':2,'x':L,'y':0,'loads':{'fy':-P}}],
        'elements':[{'id':1,'i':1,'j':2,'EI':EI,'EA':2e7,'p':3}]}
    r=solve_model(Model.model_validate(m))
    assert r['nodes'][1]['displacements']['uy']==pytest.approx(-P*L**3/(3*EI))

def test_cantilever_frame3d_tip_force():
    L=4.;EIy=2e5;P=1000.
    m={'kind':'frame3d','nodes':[{'id':1,'x':0,'y':0,'z':0,'fixed':['ux','uy','uz','rx','ry','rz']},
        {'id':2,'x':L,'y':0,'z':0,'loads':{'fz':-P}}],
        'elements':[{'id':1,'i':1,'j':2,'EIy':EIy,'EIz':5e5,'EA':2e7,'GJ':1e5,'p':3}]}
    r=solve_model(Model.model_validate(m))
    assert r['nodes'][1]['displacements']['uz']==pytest.approx(-P*L**3/(3*EIy))

def test_mechanism_detected():
    m={'kind':'truss2d','nodes':[{'id':1,'x':0,'fixed':['ux','uy']},
        {'id':2,'x':4,'loads':{'fy':-1000}}],
        'elements':[{'id':1,'i':1,'j':2,'EA':2e7,'p':1}]}
    with pytest.raises(FemError):solve_model(Model.model_validate(m))

def test_transform_symmetry():
    for kind in EXAMPLES:
        r=solve_model(Model.model_validate(EXAMPLES[kind]))
        K=np.array(r['K'])
        assert np.allclose(K,K.T,atol=1e-7)

def test_symbolic_manual_has_matrices():
    data=derive(Model.model_validate(EXAMPLES['truss2d']),'2')
    assert len(data['steps'])==10
    assert 'T' in str(data['steps'][2])
    assert any('math' in eq for step in data['steps'] for eq in step['equations'])

def test_browser_api():
    client=TestClient(app)
    assert client.get('/').status_code==200
    assert client.get('/api/examples').status_code==200
    r=client.post('/api/solve',json=EXAMPLES['beam1d']);assert r.status_code==200
    r=client.post('/api/derive',json={'model':EXAMPLES['beam1d'],'element_id':'1'});assert r.status_code==200
    assert r.json()['steps'][0]['title'].startswith('1.')

# Lagrange multiplier method: the actual solver must retain all constrained DOFs,
# construct H, and solve [K H.T; H 0][d; lambda] = [F; g].
@pytest.mark.parametrize('kind',list(EXAMPLES))
def test_lagrange_augmented_equilibrium_and_constraint_matrix(kind):
    r=solve_model(Model.model_validate(EXAMPLES[kind]))
    K=np.asarray(r['K']);F=np.asarray(r['F']);H=np.asarray(r['H']);
    g=np.asarray(r['g']);d=np.asarray(r['u']);lam=np.asarray(r['lambda'])
    n=len(d);nc=len(lam)
    assert r['solver_method']=='lagrange_multiplier'
    assert H.shape==(nc,n)
    assert nc==len(r['fixed_indices'])
    assert np.allclose(H.sum(axis=1),1)
    assert np.allclose(H[:,r['fixed_indices']],np.eye(nc))
    assert np.allclose(H@d,g,atol=1e-9)
    augmented=np.block([[K,H.T],[H,np.zeros((nc,nc))]])
    residual=augmented@np.r_[d,lam]-np.r_[F,g]
    assert np.linalg.norm(residual)<1e-5
    assert np.allclose(np.asarray(r['R']),-H.T@lam,atol=1e-5)


def test_lagrange_reaction_sign_for_indeterminate_bar():
    r=solve_model(Model.model_validate(EXAMPLES['bar1d']))
    # End support reactions are -5000 N, so with +lambda.T(Hd-g),
    # the corresponding Lagrange multipliers must be +5000 N.
    assert r['lambda']==pytest.approx([5000,5000],abs=1e-5)
    assert r['R'][r['fixed_indices'][0]]==pytest.approx(-5000)
    assert r['R'][r['fixed_indices'][1]]==pytest.approx(-5000)


def test_lagrange_manual_order_and_real_h():
    m=Model.model_validate(EXAMPLES['truss2d'])
    r=solve_model(m)
    d=derive(m,'2')
    titles=[x['title'].lower() for x in d['steps']]
    assert len(titles)==10
    assert 'element stiffness' in titles[1]
    assert 'transformation' in titles[2]
    assert 'global stiffness' in titles[3]
    assert 'boundary-condition matrix h' in titles[4]
    assert 'lagrangian' in titles[5]
    assert 'lagrange multiplier system' in titles[6]
    assert 'reactions' in titles[8]
    assert 'h' in str(d['steps'][4]['equations']).lower()
    assert d['constraint_dofs']==[r['dofs'][i] for i in r['fixed_indices']]
    assert any('<math' in eq for step in d['steps'] for eq in step['equations'])


def test_beam_reversed_element_connectivity():
    # Reversing element endpoint IDs must not change the physical response.
    base={'kind':'beam1d','nodes':[{'id':1,'x':0,'fixed':['uy','rz']},
           {'id':2,'x':4,'loads':{'fy':-1000}}],
          'elements':[{'id':1,'i':1,'j':2,'EI':2e5,'p':3}]}
    reverse=copy.deepcopy(base)
    reverse['elements'][0]['i'],reverse['elements'][0]['j']=2,1
    a=solve_model(Model.model_validate(base))
    b=solve_model(Model.model_validate(reverse))
    assert b['nodes'][1]['displacements']['uy']==pytest.approx(a['nodes'][1]['displacements']['uy'])
    assert b['nodes'][1]['displacements']['rz']==pytest.approx(a['nodes'][1]['displacements']['rz'])


def test_lagrange_system_matches_independent_reduced_solution():
    # Independent reference calculation using Kff only for verification.
    # Production code is required to use the augmented Lagrange system.
    r=solve_model(Model.model_validate(EXAMPLES['frame2d']))
    K=np.asarray(r['K']);F=np.asarray(r['F']);free=r['free_indices']
    u_ref=np.zeros(len(F))
    u_ref[free]=np.linalg.solve(K[np.ix_(free,free)],F[free])
    assert np.allclose(r['u'],u_ref,atol=1e-8)


def test_cantilever_beam_uniform_load_with_lagrange_multipliers():
    w=1200.;L=3.;EI=4e5
    m={'kind':'beam1d','nodes':[{'id':1,'x':0,'fixed':['uy','rz']},
         {'id':2,'x':L}],
       'elements':[{'id':1,'i':1,'j':2,'EI':EI,'p':3,'qy':-w}]}
    r=solve_model(Model.model_validate(m))
    assert r['nodes'][1]['displacements']['uy']==pytest.approx(-w*L**4/(8*EI))
    assert r['nodes'][1]['displacements']['rz']==pytest.approx(-w*L**3/(6*EI))
    assert r['nodes'][0]['reactions']['fy']==pytest.approx(w*L)
    assert r['nodes'][0]['reactions']['mz']==pytest.approx(w*L**2/2)


def test_indeterminate_axial_bar_mixed_p_h_and_lagrange():
    m=copy.deepcopy(EXAMPLES['bar1d'])
    m['elements'][0].update({'p':5,'h':3})
    m['elements'][1].update({'p':2,'h':2})
    r=solve_model(Model.model_validate(m))
    assert r['nodes'][1]['displacements']['ux']==pytest.approx(0.0005)
    assert r['nodes'][0]['reactions']['fx']==pytest.approx(-5000)
    assert r['nodes'][2]['reactions']['fx']==pytest.approx(-5000)
    assert len(r['elements'])==5
    assert np.linalg.norm(np.asarray(r['H'])@np.asarray(r['u'])-np.asarray(r['g']))<1e-8
