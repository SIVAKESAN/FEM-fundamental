"""Stable teaching examples, including indeterminate bars/beams/frames."""
def base(kind,nodes,elements):return {'kind':kind,'nodes':nodes,'elements':elements}

def nd(id,x,y=0,z=0,fixed=(),loads=None):
    return {'id':id,'x':x,'y':y,'z':z,'fixed':list(fixed),'loads':loads or {}}

def e(id,i,j,**kwargs):
    return {'id':id,'i':i,'j':j,'EA':2e7,'EI':2e5,'EIy':2e5,'EIz':2e5,'GJ':1e5,'p':3,'h':1,'qx':0,'qy':0,'qz':0,**kwargs}

EXAMPLES={
 'bar1d':base('bar1d',[
     nd(1,0,fixed=['ux']),nd(2,2,loads={'fx':10000}),nd(3,4,fixed=['ux'])],
     [e(1,1,2),e(2,2,3)]),
 'truss2d':base('truss2d',[
     nd(1,0,0,fixed=['ux','uy']),nd(2,4,0,fixed=['ux','uy']),nd(3,2,3,loads={'fy':-10000})],
     [e(1,1,2),e(2,1,3),e(3,2,3)]),
 'truss3d':base('truss3d',[
     nd(1,0,0,0,fixed=['ux','uy','uz']),nd(2,4,0,0,fixed=['ux','uy','uz']),
     nd(3,0,4,0,fixed=['ux','uy','uz']),nd(4,1.5,1.5,3,loads={'fz':-10000})],
     [e(1,1,4),e(2,2,4),e(3,3,4)]),
 'beam1d':base('beam1d',[
     nd(1,0,fixed=['uy','rz']),nd(2,2),nd(3,4,fixed=['uy','rz'])],
     [e(1,1,2,qy=-2000),e(2,2,3,qy=-2000)]),
 'frame2d':base('frame2d',[
     nd(1,0,0,fixed=['ux','uy','rz']),nd(2,0,3),nd(3,4,3,loads={'fy':-5000}),
     nd(4,4,0,fixed=['ux','uy','rz'])],
     [e(1,1,2),e(2,2,3,qy=-1000),e(3,3,4)]),
 'frame3d':base('frame3d',[
     nd(1,0,0,0,fixed=['ux','uy','uz','rx','ry','rz']),
     nd(2,0,0,3,loads={'fx':2000,'fy':-1000,'fz':-5000}),
     nd(3,4,0,3),
     nd(4,4,0,0,fixed=['ux','uy','uz','rx','ry','rz'])],
     [e(1,1,2),e(2,2,3),e(3,3,4)])
}
