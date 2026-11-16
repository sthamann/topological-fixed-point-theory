// Exact audit of the E8 quarter shift in the native colour frame.
// Coordinates are doubled for roots and in twelfths for torus shifts.
const basis = [
  [1,-1,-1,-1,-1,-1,-1,1], [2,2,0,0,0,0,0,0], [-2,2,0,0,0,0,0,0],
  [0,-2,2,0,0,0,0,0], [0,0,-2,2,0,0,0,0], [0,0,0,-2,2,0,0,0],
  [0,0,0,0,-2,2,0,0], [0,0,0,0,0,-2,2,0]
];
const colourRoots = [[0,0,0,0,0,2,-2,0], [0,0,0,0,0,0,2,-2]];
const quarter = [0,0,0,0,0,6,6,6];
const geometricQuarter = [0,0,3,3,6,6,9,9];
const colourCentre = [0,0,0,0,0,8,-4,-4];
const nativeCounts = [60,64,60,64];
const zero = () => Array(8).fill(0);
const scale = (n,v) => v.map(x => n*x);
const add = (v,w) => v.map((x,i) => x+w[i]);
const dot = (v,w) => v.reduce((s,x,i) => s+x*w[i],0);
const mod = (a,n) => ((a%n)+n)%n;
const check = (ok,why) => { if (!ok) throw Error(why); };

const roots = [];
for (let i=0; i<8; i++) for (let j=i+1; j<8; j++) {
  for (const x of [-2,2]) for (const y of [-2,2]) {
    const v=zero();
    v[i]=x;
    v[j]=y;
    roots.push(v);
  }
}
for (let mask=0; mask<256; mask++) {
  const v=zero().map((_,i)=>(mask>>i&1) ? -1 : 1);
  if (v.filter(x=>x===-1).length%2===0) roots.push(v);
}
check(roots.length===240 && new Set(roots.map(String)).size===240, "E8 root census");
const gram=basis.map(v=>basis.map(w=>dot(v,w)/4));
// The native simple roots form the full lattice, not a sublattice of odd index.
const elimination=gram.map(row=>row.slice());
let prior=1;
for(let k=0;k<7;k++) {
  const pivot=elimination[k][k];
  check(pivot>0,"positive leading E8 Gram minor");
  for(let i=k+1;i<8;i++) for(let j=k+1;j<8;j++)
    elimination[i][j]=(pivot*elimination[i][j]-elimination[i][k]*elimination[k][j])/prior;
  prior=pivot;
}
check(elimination[7][7]===1,"unimodular E8 basis");

function phaseSearch() {
  const compatible=Array(4).fill(0);
  const matches=Array(4).fill(0);
  const typeZeroSpectra={};
  let fixesColour=0;
  for(let id=0; id<65536; id++) {
    let v=zero();
    for (let i=0;i<8;i++) v=add(v,scale(id>>(2*i)&3,basis[i]));
    const square=dot(v,v);
    check(square%8===0,"even E8 lattice");
    const type=square/8%4;
    const counts=[8,0,0,0];
    for(const r of roots) {
      const pairing=dot(v,r);
      check(pairing%4===0,"integral root pairing");
      counts[mod(pairing/4,4)]++;
    }
    const fixed=colourRoots.every(r=>mod(dot(r,v),16)===0);
    if(fixed) {
      fixesColour++;
      if(type===0 && counts[1]>0) typeZeroSpectra[counts]=(typeZeroSpectra[counts]??0)+1;
    }
    if(counts.every((n,i)=>n===nativeCounts[i])) {
      matches[type]++;
      if(fixed) compatible[type]++;
    }
  }
  check(String(matches)==="0,0,15120,0" && String(compatible)==="0,0,594,0" &&
    fixesColour===4096 && typeZeroSpectra["64,64,56,64"]===864 &&
    Object.keys(typeZeroSpectra).length===1,"complete clock spectrum census");
  return {enumerated:65536, preservesColourRoots:fixesColour,
    matchingNativeByType:matches, matchingNativeAndColourByType:compatible, typeZeroSpectra};
}

function belongsToE8(v) {
  const residue=mod(v[0],12);
  return (residue===0 || residue===6) && v.every(x=>mod(x,12)===residue) &&
    mod(v.reduce((a,b)=>a+b,0),24)===0;
}

function shortestPoints(shift) {
  const sectors=[];
  for (const offset of [0,6]) {
    let states=[{norm:0,points:[[]]},null];
    for (const x of shift) {
      const choices=[0,1].map(parity=>{
        const start=x+offset+12*parity;
        const k=Math.floor(-start/24);
        const values=[start+24*k,start+24*(k+1)];
        const square=Math.min(...values.map(v=>v*v));
        return {square,values:values.filter(v=>v*v===square)};
      });
      const next=[null,null];
      for(let parity=0;parity<2;parity++) {
        if(!states[parity]) continue;
        for(let step=0;step<2;step++) {
          const end=parity^step;
          const squared=states[parity].norm+choices[step].square;
          if(!next[end] || squared<next[end].norm) next[end]={norm:squared,points:[]};
          if(squared===next[end].norm) {
            for(const p of states[parity].points)
              for(const value of choices[step].values) next[end].points.push([...p,value]);
          }
        }
      }
      states=next;
    }
    sectors.push(states[0]);
  }
  const norm=Math.min(...sectors.map(s=>s.norm));
  return {norm,points:sectors.filter(s=>s.norm===norm).flatMap(s=>s.points)};
}

function cyclicType(n,shift) {
  check(belongsToE8(scale(n,shift)),"declared order does not annihilate the shift");
  for(let i=1;i<n;i++) check(!belongsToE8(scale(i,shift)),"smaller shift order");
  const scaled=n*n*shortestPoints(shift).norm/288;
  check(Number.isInteger(scaled),"nonintegral cyclic type");
  return scaled%n;
}

function twistedSectors() {
  const sectors=[0,1,2,3].map(power=>{
    const {norm,points}=shortestPoints(scale(power,quarter));
    const weights=points.map(p=>colourRoots.map(r=>dot(p,r)/24));
    return {power, conformalWeight:norm/288, multiplicity:points.length,weights,points};
  });
  check(sectors.map(s=>s.conformalWeight).join()==="0,0.375,0.5,0.375" &&
        sectors.map(s=>s.multiplicity).join()==="1,4,16,4","twist energies and degeneracies");
  check(sectors[1].weights.map(String).sort().join(";")==="-1,1;0,-1;0,0;1,0","triplet plus singlet");
  check(sectors[3].weights.map(w=>String(scale(-1,w))).sort().join(";")==="-1,1;0,-1;0,0;1,0","conjugate triplet");
  check(sectors[2].points.every(p=>p.reduce((n,x)=>n+Math.abs(x),0)===12 &&
    p.every(x=>x%12===0)),"half turn vector charges");
  for(const sector of [sectors[1],sectors[3]]) {
    const present=new Set(sector.points.map(String));
    for(const r of colourRoots) {
      check(sector.points.filter(p=>present.has(String(add(p,scale(6,r))))).length===1,
        "nonzero colour step on the twisted triplet");
    }
  }
  return sectors.map(({points,...report})=>report);
}

function backgroundResponse() {
  const clock=scale(4/12,quarter);
  const colour=scale(3/12,colourCentre);
  check(dot(clock,clock)===12 && dot(clock,colour)===0 && dot(colour,colour)===6,
    "clock and colour inflow matrix");
  return {clockCharge:basis.map(r=>dot(r,clock)/2),colourCharge:basis.map(r=>dot(r,colour)/2),
    responseUpToOrientation:[[12,0],[0,6]], discreteTypes:[2,0]};
}

function discriminantGlue() {
  const q8=(a,b)=>mod(3*a*a+5*b*b,8);
  const subgroups=new Map();
  for(let a=0;a<4;a++) for(let b=0;b<4;b++) {
    const group=[0,1,2,3].map(k=>[mod(k*a,4),mod(k*b,4)]);
    if(new Set(group.map(String)).size===4 && group.every(p=>q8(...p)===0)) {
      subgroups.set(group.map(String).sort().join(";"),group);
    }
  }
  const lagrangians=[...subgroups.values()];
  check(lagrangians.length===2 && q8(2,0)===4 && q8(0,2)===4,
    "D5+A3 isotropic subgroups");
  const shared=lagrangians[0].filter(p=>lagrangians[1].some(q=>String(p)===String(q)));
  check(shared.map(String).sort().join(";")==="0,0;2,2",
    "no transverse isotropic quarter extension");
  return {lagrangians,shared,transverseIsotropicComplement:false};
}

function audit() {
  const types=[cyclicType(4,quarter),cyclicType(4,geometricQuarter),cyclicType(2,scale(2,quarter)),
    cyclicType(3,colourCentre),cyclicType(12,add(quarter,colourCentre))];
  check(types.join()==="2,2,0,0,6","cyclic types");
  const omegaExponent=(a,b,c)=>types[0]*a*Math.floor((b+c)/4);
  let cocycleExponent=0;
  for(let j=0;j<4;j++) cocycleExponent+=omegaExponent(1,j,1);
  check(cocycleExponent%4===2,"nontrivial cyclic cocycle invariant");
  for(let a=0;a<4;a++) for(let b=0;b<4;b++) for(let c=0;c<4;c++) for(let d=0;d<4;d++)
    check(mod(omegaExponent(b,c,d)-omegaExponent((a+b)%4,c,d)+omegaExponent(a,(b+c)%4,d)
      -omegaExponent(a,b,(c+d)%4)+omegaExponent(a,b,c),4)===0,"cyclic cocycle equation");
  const vector=[12,0,0,0,0,0,0,0], spinor=Array(8).fill(6);
  check(belongsToE8(add(scale(2,quarter),vector)) && belongsToE8(spinor) &&
    mod(dot(vector,spinor),144)===72,"vector and spinor half-integer pairing");
  return {types:{quarter:types[0],geometricQuarter:types[1],halfTurn:types[2],
    colourCentre:types[3],joint:types[4]}, quarticCocycleProduct:-1,
    phaseSearch:phaseSearch(),groundSectors:twistedSectors(),vectorSpinorMonodromy:-1,
    inflow:backgroundResponse(), discriminantGlue:discriminantGlue(), physicalR113SourceDerived:false};
}

console.log(JSON.stringify(audit(),null,2));
