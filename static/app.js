/* FEM Learning Laboratory — no npm or JavaScript frameworks needed. */
'use strict';
const $ = (id)=>document.getElementById(id);
const typeNames={bar1d:'1D BAR (axial)',truss2d:'2D TRUSS',truss3d:'3D TRUSS',beam1d:'1D BEAM (bending)',frame2d:'2D FRAME',frame3d:'3D SPACE FRAME'};
const typedofs={bar1d:['ux'],truss2d:['ux','uy'],truss3d:['ux','uy','uz'],beam1d:['uy','rz'],frame2d:['ux','uy','rz'],frame3d:['ux','uy','uz','rx','ry','rz']};
const forceKey={ux:'fx',uy:'fy',uz:'fz',rx:'mx',ry:'my',rz:'mz'};
const minOrder=(kind)=>['beam1d','frame2d','frame3d'].includes(kind)?3:1;
let examples={},model=null,result=null,selectedElement=null,showManual=false;
const copy=(o)=>JSON.parse(JSON.stringify(o));
const escapeHtml=(v)=>String(v??'').replace(/[&<>"']/g,s=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[s]));
const num=(n,d=5)=>{let v=Number(n);if(!Number.isFinite(v))return '—';if(v===0)return '0';return (Math.abs(v)<.0001||Math.abs(v)>1e6)?v.toExponential(3):v.toLocaleString('en-US',{maximumFractionDigits:d})};
const numberInput=(value,key,index,group,width='')=>`<input aria-label="${key} ${index+1}" type="number" step="any" class="${width}" data-group="${group}" data-index="${index}" data-field="${key}" value="${escapeHtml(value??0)}">`;
function setStatus(text,error=false,success=false){$('status').textContent=text;$('status').className='instructions'+(error?' error':success?' success':'');}
function invalidate(){result=null;showManual=false;$('manualButton').disabled=true;$('downloadCode').disabled=true;$('downloadMatrices').disabled=true;$('print').disabled=true;$('solutionSection').hidden=true;$('overview').className='overview empty';$('overview').innerHTML='Model changed. Press <strong>SOLVE STRUCTURE</strong> to recalculate.';$('tables').hidden=true;draw();}
function supportPresets(){let ds=typedofs[model.kind];let pin=ds.filter(x=>x[0]==='u');let roller=ds.filter(x=>x==='uy');if(model.kind==='bar1d')roller=['ux'];if(model.kind==='beam1d')pin=['uy'];return {Free:[],Fixed:ds,Pin:pin,Roller:roller};}
function supportValue(node){const fixed=[...node.fixed].sort().join(',');const names=supportPresets();return Object.entries(names).find(([k,v])=>[...v].sort().join(',')===fixed)?.[0]||'Custom'}
function resetForType(kind){model=copy(examples[kind]);result=null;selectedElement=String(model.elements[0]?.id||'1');render();invalidate();setStatus(`Example loaded: ${typeNames[kind]}. Edit numerical inputs, or press SOLVE STRUCTURE.`);}
function render(){if(!model)return;$('kind').value=model.kind;renderNodes();renderElements();draw();}
function renderNodes(){let dofs=typedofs[model.kind];let fs=dofs.map(x=>forceKey[x]);let fields=['x',...(model.kind==='bar1d'||model.kind==='beam1d'?[]:['y']),...(model.kind==='truss3d'||model.kind==='frame3d'?['z']:[])];
  $('nodesTable').innerHTML=`<thead><tr><th>ID</th>${fields.map(x=>`<th>${x} (m)</th>`).join('')}<th>Support</th><th>Restrained DOFs</th>${fs.map(x=>`<th>${x}</th>`).join('')}<th></th></tr></thead><tbody>`+
    model.nodes.map((n,i)=>`<tr><td>N${n.id}</td>${fields.map(k=>`<td>${numberInput(n[k],k,i,'node', 'short')}</td>`).join('')}
      <td><select data-action="support" data-index="${i}">${[...Object.keys(supportPresets()),'Custom'].map(v=>`<option ${supportValue(n)===v?'selected':''}>${v}</option>`).join('')}</select></td>
      <td><input aria-label="restrained DOFs node ${n.id}" class="wide" data-group="node" data-index="${i}" data-field="fixed" type="text" value="${escapeHtml(n.fixed.join(', '))}" placeholder="ux,uy"></td>
      ${fs.map(key=>`<td>${numberInput(n.loads[key]||0,key,i,'load','short')}</td>`).join('')}<td><button class="delete" title="Delete node" data-action="delete-node" data-index="${i}">×</button></td></tr>`).join('')+'</tbody>';
}
function renderElements(){let kind=model.kind;let props=['EA',...(['beam1d','frame2d'].includes(kind)?['EI']:[]),...(kind==='frame3d'?['EIy','EIz','GJ']:[])];let qf=['beam1d','frame2d'].includes(kind)?['qy']:kind==='frame3d'?['qx','qy','qz']:[];
  $('elementsTable').innerHTML=`<thead><tr><th>ID</th><th>i</th><th>j</th>${props.map(x=>`<th>${x}</th>`).join('')}<th>p</th><th>h</th>${qf.map(x=>`<th>${x}</th>`).join('')}<th></th></tr></thead><tbody>`+
    model.elements.map((e,i)=>`<tr><td>E${e.id}</td>${['i','j'].map(k=>`<td>${numberInput(e[k],k,i,'element','short')}</td>`).join('')}
      ${props.map(k=>`<td>${numberInput(e[k],k,i,'element')}</td>`).join('')}
      <td><input title="Polynomial order" type="number" min="${minOrder(kind)}" max="6" step="1" class="short" data-group="element" data-index="${i}" data-field="p" value="${e.p}"></td>
      <td><input title="h subdivisions" type="number" min="1" max="8" step="1" class="short" data-group="element" data-index="${i}" data-field="h" value="${e.h}"></td>
      ${qf.map(k=>`<td>${numberInput(e[k]||0,k,i,'element')}</td>`).join('')}
      <td><button class="delete" data-action="delete-element" data-index="${i}" title="Delete element">×</button></td></tr>`).join('')+'</tbody>';
}
function updateField(el){let index=+el.dataset.index,field=el.dataset.field,group=el.dataset.group;
  const object=group==='element'?model.elements[index]:model.nodes[index];if(!object)return;
  if(field==='fixed')object.fixed=el.value.split(',').map(v=>v.trim().toLowerCase()).filter(Boolean);
  else if(group==='load'){object.loads[field]=Number(el.value)||0;}
  else object[field]=Number(el.value)||0;
  invalidate();if(field==='fixed')renderNodes();
}
function project(a){return [a.x+a.z*.42, a.y+a.z*.32];}
function draw(){if(!model)return;const W=760,H=370;
  let verts=model.nodes.map(n=>{let [x,y]=project(n);return{id:n.id,x,y,n}}),xx=verts.map(n=>n.x),yy=verts.map(n=>n.y);
  let xmin=Math.min(...xx),xmax=Math.max(...xx),ymin=Math.min(...yy),ymax=Math.max(...yy);let span=Math.max(xmax-xmin,ymax-ymin,1);
  const margin=76; const scale=Math.min((W-2*margin)/Math.max(xmax-xmin,span*.3),(H-2*margin)/Math.max(ymax-ymin,span*.3));
  const xMid=(xmin+xmax)/2,yMid=(ymin+ymax)/2;
  let pt=(v)=>[W/2+(v.x-xMid)*scale,H/2-(v.y-yMid)*scale];
  let m=new Map(verts.map(v=>[v.id,pt(v)]));
  let out=`<rect x="0" y="0" width="${W}" height="${H}" fill="white"/>`;
  for(let i=1;i<10;i++){let x=i*W/10,y=i*H/10;out+=`<line class="svg-grid" x1="${x}" y1="0" x2="${x}" y2="${H}"/><line class="svg-grid" x1="0" y1="${y}" x2="${W}" y2="${y}"/>`}
  const forceVals=model.nodes.flatMap(n=>Object.values(n.loads||{})).map(Math.abs);let maxF=Math.max(...forceVals,1);
  let displacementScale=0;
  let nd=new Map();
  if(result){const maxU=result.summary.max_translation_m||0;displacementScale=maxU>0?Math.min(span*.22/maxU,2e7):0; for(let n of result.nodes)nd.set(n.id,n.displacements);}
  const original=new Map(model.nodes.map(n=>[n.id,n]));
  for(let e of model.elements){let a=m.get(e.i),b=m.get(e.j);if(!a||!b)continue;
    let selected=String(e.id)===String(selectedElement)?' selected':'';
    out+=`<line class="svg-member${selected}" data-element="${e.id}" x1="${a[0]}" y1="${a[1]}" x2="${b[0]}" y2="${b[1]}"/>`;
    for(let j=1;j<e.h;j++){let r=j/e.h;let px=a[0]+(b[0]-a[0])*r,py=a[1]+(b[1]-a[1])*r;
      out+=`<circle cx="${px}" cy="${py}" r="3.1" fill="white" stroke="black" stroke-width="1.4" pointer-events="none"/>`;
    }
    const mid=[(a[0]+b[0])/2,(a[1]+b[1])/2];out+=`<text class="element-id" x="${mid[0]+4}" y="${mid[1]-9}">E${e.id} / p${e.p} h${e.h}</text>`;
    if(result&&$('showDeformed').checked&&nd.has(e.i)&&nd.has(e.j)){
      let n1=original.get(e.i),n2=original.get(e.j);let t1=nd.get(e.i),t2=nd.get(e.j);
      let sh=(n,d)=>{const x=n.x+(d.ux||0)*displacementScale,y=n.y+(d.uy||0)*displacementScale,z=n.z+(d.uz||0)*displacementScale;return pt({x:x+z*.42,y:y+z*.32})};
      let p1=sh(n1,t1),p2=sh(n2,t2);
      out+=`<line class="svg-deformed" x1="${p1[0]}" y1="${p1[1]}" x2="${p2[0]}" y2="${p2[1]}"/>`;
    }
  }
  for(let node of model.nodes){let a=m.get(node.id);if(!a)continue;let [x,y]=a;
    if(node.fixed?.length){const fully=node.fixed.length===typedofs[model.kind].length;
      if(fully){out+=`<line class="svg-support" x1="${x-15}" y1="${y+10}" x2="${x+15}" y2="${y+10}"/>`;
      for(let v=-12;v<=12;v+=6)out+=`<line class="svg-hatch" x1="${x+v}" y1="${y+10}" x2="${x+v-7}" y2="${y+18}"/>`;
      }else{out+=`<path class="svg-support" d="M${x} ${y+5} l-14 16 h28 Z M${x-18} ${y+24} h36"/>`}
    }
    out+=`<circle class="svg-node" cx="${x}" cy="${y}" r="5"/>`;
    if($('showNodes').checked)out+=`<text class="svg-label" x="${x+8}" y="${y+19}">${node.id}</text>`;
    // nodal forces in X/Y shown as restrained monochrome arrows
    let fx=node.loads?.fx||0,fy=node.loads?.fy||0,fz=node.loads?.fz||0;
    function arrow(dx,dy,value,label){let length=28+Math.min(13,Math.abs(value)/maxF*13);let s=value>=0?1:-1;let vx=dx*s,vy=dy*s;
      let ex=x+vx*length,ey=y+vy*length;
      out+=`<path class="svg-load" d="M${x+vx*8} ${y+vy*8} L${ex} ${ey} m${-vx*8-vy*3} ${-vy*8+vx*3} l${vx*8+vy*3} ${vy*8-vx*3} l${-vx*8+vy*3} ${-vy*8-vx*3}"/>`;
      out+=`<text class="svg-loadtext" x="${ex+vx*10+4}" y="${ey+vy*10}">${label} ${num(value,2)}</text>`;
    }
    if(fx)arrow(1,0,fx,'Fx');if(fy)arrow(0,-1,fy,'Fy');if(fz)arrow(.7,-.7,fz,'Fz');
  }
  out+=`<text class="svg-loadtext" x="15" y="22">${model.kind.toUpperCase()} • SI UNITS</text>`;
  out+=`<path class="svg-load" d="M50 335 H95 M90 330 L95 335 L90 340 M50 335 V294 M45 299 L50 294 L55 299"/><text class="svg-loadtext" x="97" y="338">x</text><text class="svg-loadtext" x="41" y="285">y</text>`;
  if(model.kind.includes('3d'))out+=`<text class="svg-loadtext" x="16" y="355">Oblique projection of 3D coordinates</text>`;
  $('drawing').innerHTML=out;
}
function renderResults(){let s=result.summary;let metrics=[['MODEL',typeNames[result.kind]],['TOTAL DEGREES OF FREEDOM',num(result.n_total_dofs)],['ELEMENTS AFTER h-REFINEMENT',num(s.elements_after_h_refinement)],['MAX. TRANSLATION',num(s.max_translation_m*1e3)+' mm'],['MAX. ROTATION',num(s.max_rotation_rad,6)+' rad'],['FREE-DOF RESIDUAL',num(s.free_equilibrium_residual_N,6)],['SOLVER','Lagrange multipliers'],['CONSTRAINT ERROR',num(s.constraint_residual,6)],['CONSTRAINTS',num(result.fixed_indices.length)]];
  $('overview').className='overview';$('overview').innerHTML=metrics.map(([a,b])=>`<div class="metric"><div class="tag">${escapeHtml(a)}</div><strong>${escapeHtml(b)}</strong></div>`).join('');
  let dofs=typedofs[model.kind];let forces=dofs.map(x=>forceKey[x]);
  let nodetable=`<h3>Nodal displacements and support reactions</h3><table class="output-table"><thead><tr><th>Node</th>${dofs.map(k=>`<th>${k} (${k.startsWith('r')?'rad':'m'})</th>`).join('')}${forces.map(k=>`<th>R${k.slice(1)} (${k.startsWith('m')?'N·m':'N'})</th>`).join('')}</tr></thead><tbody>`+
  result.nodes.map(n=>`<tr><td>${n.generated?'*':''}${n.id}</td>${dofs.map(d=>`<td>${num(n.displacements[d],7)}</td>`).join('')}${forces.map(f=>`<td>${n.reactions[f]===undefined?'—':num(n.reactions[f])}</td>`).join('')}</tr>`).join('')+'</tbody></table>';
  let elemResults=`<h3>Element end-force vectors (local axes, SI units)</h3><table class="output-table"><thead><tr><th>Element</th><th>Nodes</th><th>Length (m)</th><th>End forces (local DOF order)</th></tr></thead><tbody>`+
  result.elements.map(e=>`<tr><td>${escapeHtml(e.id)}</td><td>${e.i}–${e.j}</td><td>${num(e.L)}</td><td>${escapeHtml(e.end_forces_local.map(x=>num(x,4)).join(', '))}</td></tr>`).join('')+'</tbody></table>';
  let multiplierTable=`<h3>Lagrange multiplier boundary conditions (H·d = g)</h3><table class="output-table"><thead><tr><th>Constraint</th><th>Restrained DOF</th><th>λ</th><th>Reaction = −λ</th></tr></thead><tbody>`+
    result.fixed_indices.map((idx,i)=>`<tr><td>${i+1}</td><td>${escapeHtml(result.dofs[idx])}</td><td>${num(result.lambda[i],7)}</td><td>${num(-result.lambda[i],7)}</td></tr>`).join('')+'</tbody></table>';
  $('tables').innerHTML=`<div class="tables-wrap">${nodetable}<p class="hint">* Generated h-refinement node. Reaction entries appear only at restrained DOFs.</p>${multiplierTable}${elemResults}<h3>Global degree-of-freedom sequence</h3><p class="hint">${escapeHtml(result.dofs.map((s,i)=>`${i+1}:${s}`).join(' · '))}</p></div>`;
  $('tables').hidden=false;$('manualButton').disabled=false;$('downloadCode').disabled=false;$('downloadMatrices').disabled=false;$('print').disabled=false;draw();
}
async function solve(){if(!model)return;setStatus('Solving FEM equations using Python backend…');$('solve').disabled=true;
 try{const response=await fetch('/api/solve',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(model)});
  const payload=await response.json();if(!response.ok)throw Error(typeof payload.detail==='string'?payload.detail:JSON.stringify(payload.detail));
  result=payload;renderResults();setStatus(`Solved by Lagrange multipliers: ${result.n_total_dofs} DOFs, ${result.fixed_indices.length} constraints. Click SHOW COMPLETE SOLUTION for Kₑ, Kᴳ, H and the augmented system.`,false,true);
 }catch(err){setStatus('Analysis error: '+err.message,true);result=null;$('manualButton').disabled=true;draw();}
 finally{$('solve').disabled=false;}}
function selectMemberById(id){selectedElement=String(id);draw();if(showManual){$('selectedElement').value=String(id);renderManual();}}
function buildSelect(){let all=model.elements.map(e=>[String(e.id),`Element ${e.id} (${e.i} → ${e.j})`]);
  $('selectedElement').innerHTML=all.map(([v,t])=>`<option value="${escapeHtml(v)}" ${v===String(selectedElement)?'selected':''}>${escapeHtml(t)}</option>`).join('');
  if(!all.some(a=>a[0]===String(selectedElement)))selectedElement=String(all[0]?.[0]||'1');$('selectedElement').value=selectedElement;}
async function showSolution(){if(!result){setStatus('Solve the structure before opening the solution.',true);return;}showManual=true;$('solutionSection').hidden=false;buildSelect();await renderManual();$('solutionSection').scrollIntoView({behavior:'smooth',block:'start'});}
async function renderManual(){if(!model||!result)return;let index=selectedElement;let box=$('derivation');box.innerHTML='<p>Computing the symbolic step-by-step derivation in Python…</p>';
 try{let response=await fetch('/api/derive',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({model,element_id:String(index)})});let data=await response.json();if(!response.ok)throw Error(data.detail||'Derivation request failed');
    box.innerHTML=`<div class="manual-note"><strong>Selected member:</strong> ${escapeHtml(data.element)} &nbsp;|&nbsp; ${escapeHtml(typeNames[data.kind])} &nbsp;|&nbsp; Length = ${num(data.selected_geometry.length_m)} m &nbsp;|&nbsp; p = ${data.selected_geometry.p} &nbsp;|&nbsp; h = ${data.selected_geometry.h}<br/><strong>Element DOFs:</strong> ${escapeHtml(data.element_dofs.join(', '))}<br/><strong>Solution method:</strong> Lagrange multipliers · Kₑ → Kᴳ → H → augmented system → d, λ → reactions</div>`+
      data.steps.map((s,i)=>`<details ${i<=2?'open':''}><summary>${escapeHtml(s.title)}</summary><div class="step-body"><p>${escapeHtml(s.body)}</p>${s.equations.map(eq=>`<div class="math-container">${eq}</div>`).join('')}${s.note?`<div class="manual-note">${escapeHtml(s.note)}</div>`:''}</div></details>`).join('');
  }catch(e){box.innerHTML=`<p class="instructions error">${escapeHtml(e.message)}</p>`;}}
function addNode(){if(!model)return;let id=Math.max(...model.nodes.map(n=>n.id),0)+1;let last=model.nodes[model.nodes.length-1];model.nodes.push({id,x:(last?.x||0)+2,y:last?.y||0,z:last?.z||0,fixed:[],loads:{}});render();invalidate();}
function addElement(){if(!model||model.nodes.length<2)return;let id=Math.max(...model.elements.map(e=>e.id),0)+1;let n=model.nodes;
  model.elements.push({id,i:n[n.length-2].id,j:n[n.length-1].id,EA:2e7,EI:2e5,EIy:2e5,EIz:2e5,GJ:1e5,p:Math.max(minOrder(model.kind),3),h:1,qx:0,qy:0,qz:0});render();invalidate();}
function download(name,content,type){let url=URL.createObjectURL(new Blob([content],{type}));let a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
function saveModel(){download('fem-model.json',JSON.stringify(model,null,2),'application/json');}
function downloadMatrices(){if(!result)return;const payload={
  kind:result.kind,solver_method:result.solver_method,dof_order:result.dofs,
  element_stiffness_note:'For each selected member, inspect Ke and T in SHOW COMPLETE SOLUTION.',
  K_G:result.K,F:result.F,H:result.H,g:result.g,
  augmented_matrix:'A = [[K_G, transpose(H)], [H, zeros]]',
  d:result.u,lambda:result.lambda,reactions:result.R,
  constraint_indices_zero_based:result.fixed_indices,
  reaction_sign_convention:'R = K_G*d - F = -transpose(H)*lambda'
};download('fem-lagrange-matrices.json',JSON.stringify(payload,null,2),'application/json');}
function downloadPython(){let content=`"""Reproducible Python FEM calculation. Run from FEM_Learning_Laboratory project root."""\nimport json\nfrom backend.models import Model\nfrom backend.solver import solve_model\n\nraw = json.loads(r'''${JSON.stringify(model,null,2)}''')\nmodel = Model.model_validate(raw)\nresult = solve_model(model)\nprint('SUMMARY:', result['summary'])\nfor node in result['nodes']:\n    print('NODE', node['id'], 'DISPLACEMENTS', node['displacements'], 'REACTIONS', node['reactions'])\nprint('GLOBAL STIFFNESS MATRIX K_G:')\nfor row in result['K']:\n    print(row)\nprint('BOUNDARY-CONDITION MATRIX H:')\nfor row in result['H']:\n    print(row)\nprint('LAGRANGE MULTIPLIERS:', result['lambda'])\nprint('GLOBAL DISPLACEMENTS:', result['u'])\n`;download('fem-solution.py',content,'text/x-python');}
async function start(){try{let r=await fetch('/api/examples');if(!r.ok)throw Error('Unable to load examples');examples=await r.json();
  $('kind').innerHTML=Object.entries(typeNames).map(([v,k])=>`<option value="${v}">${k}</option>`).join('');resetForType('bar1d');}
  catch(e){setStatus('Cannot connect to Python backend: '+e.message,true);}}
$('kind').addEventListener('change',e=>resetForType(e.target.value));
$('example').addEventListener('click',()=>resetForType(model.kind));
$('nodesTable').addEventListener('change',e=>{let t=e.target;if(t.dataset.action==='support'){let n=model.nodes[+t.dataset.index];if(t.value!=='Custom'){n.fixed=[...supportPresets()[t.value]];invalidate();renderNodes()}return;}if(t.dataset.group)updateField(t)});
$('elementsTable').addEventListener('change',e=>{if(e.target.dataset.group)updateField(e.target)});
function handleTableButton(e,which){let b=e.target.closest('button[data-action]');if(!b)return;let idx=+b.dataset.index;if(b.dataset.action==='delete-node'){
  if(model.nodes.length<=2){setStatus('At least 2 nodes are required.',true);return;}let id=model.nodes[idx].id;model.nodes.splice(idx,1);model.elements=model.elements.filter(el=>el.i!==id&&el.j!==id);if(!model.elements.length)setStatus('Add an element connecting the remaining nodes.',true);
 }else if(b.dataset.action==='delete-element'){if(model.elements.length<=1){setStatus('At least 1 element is required.',true);return;}model.elements.splice(idx,1);}render();invalidate()}
$('nodesTable').addEventListener('click',handleTableButton);
$('elementsTable').addEventListener('click',handleTableButton);
$('addNode').addEventListener('click',addNode);$('addElement').addEventListener('click',addElement);
for(let [id,key,delta,lo,hi] of [['incH','h',1,1,8],['decH','h',-1,1,8],['incP','p',1,1,6],['decP','p',-1,1,6]]){
  $(id).addEventListener('click',()=>{for(let el of model.elements){el[key]=Math.min(hi,Math.max(key==='p'?minOrder(model.kind):lo,(el[key]||lo)+delta))}render();invalidate()});}
$('solve').addEventListener('click',solve);$('manualButton').addEventListener('click',showSolution);
$('refreshDerivation').addEventListener('click',renderManual);
$('selectedElement').addEventListener('change',e=>selectMemberById(e.target.value));
$('hideManual').addEventListener('click',()=>{showManual=false;$('solutionSection').hidden=true;});
$('save').addEventListener('click',saveModel);$('load').addEventListener('click',()=>$('file').click());
$('file').addEventListener('change',async e=>{let file=e.target.files?.[0];if(!file)return;try{let data=JSON.parse(await file.text());if(!typeNames[data.kind]||!Array.isArray(data.nodes)||!Array.isArray(data.elements))throw Error('Unrecognized model format');model=data;selectedElement=String(model.elements[0]?.id||'1');render();invalidate();setStatus('Model opened. Press SOLVE STRUCTURE.');}catch(ex){setStatus('Cannot open JSON: '+ex.message,true)}e.target.value='';});
$('downloadCode').addEventListener('click',downloadPython);$('downloadMatrices').addEventListener('click',downloadMatrices);$('print').addEventListener('click',()=>window.print());
$('showNodes').addEventListener('change',draw);$('showDeformed').addEventListener('change',draw);
$('drawing').addEventListener('click',e=>{let el=e.target.closest('[data-element]');if(el)selectMemberById(el.dataset.element)});
window.addEventListener('load',start);
