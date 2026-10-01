/* Oblicua · detectar las esquinas de una superficie tocándola.
   1. Bordes del fotograma (Sobel) y líneas rectas largas (transformada de Hough, votando solo cerca de la
      orientación del borde).
   2. Con pares de líneas casi paralelas a cada lado del punto tocado se arman cuadriláteros que lo contienen.
   3. Cada uno se puntúa por cuánto borde real hay a lo largo de sus 4 lados; los mejores se afinan ajustando
      cada lado a los puntos de borde cercanos.
   Coordenadas: píxeles de la imagen de trabajo (centro del píxel i en i). */

function sobel(g,w,h){
  const gx=new Float32Array(w*h),gy=new Float32Array(w*h),mag=new Float32Array(w*h);
  for(let y=1;y<h-1;y++)for(let x=1;x<w-1;x++){
    const o=y*w+x,a=g[o-w-1],b=g[o-w],c=g[o-w+1],d=g[o-1],f=g[o+1],k=g[o+w-1],l=g[o+w],m=g[o+w+1];
    const X=(c+2*f+m)-(a+2*d+k),Y=(k+2*l+m)-(a+2*b+c);gx[o]=X;gy[o]=Y;mag[o]=Math.hypot(X,Y);
  }
  return{gx,gy,mag};
}
function samp(I,w,h,x,y){
  if(x<0)x=0;else if(x>w-1.001)x=w-1.001;if(y<0)y=0;else if(y>h-1.001)y=h-1.001;
  const ix=x|0,iy=y|0,fx=x-ix,fy=y-iy,o=iy*w+ix,a=I[o]+(I[o+1]-I[o])*fx,b=I[o+w]+(I[o+w+1]-I[o+w])*fx;return a+(b-a)*fy;
}
/** Líneas rectas: [{th,rho,votes}] con x·cos(th)+y·sin(th)=rho. */
function hough(E,w,h,thr){
  const NT=180,diag=Math.ceil(Math.hypot(w,h)),NR=2*diag+1,acc=new Uint16Array(NT*NR);
  const cs=new Float32Array(NT),sn=new Float32Array(NT);for(let t=0;t<NT;t++){cs[t]=Math.cos(t*Math.PI/NT);sn[t]=Math.sin(t*Math.PI/NT);}
  for(let y=2;y<h-2;y++)for(let x=2;x<w-2;x++){
    const o=y*w+x;if(E.mag[o]<thr)continue;
    // solo cerca de la orientación del borde (la normal de la línea es el gradiente)
    let t0=Math.round(((Math.atan2(E.gy[o],E.gx[o])+Math.PI)%Math.PI)/Math.PI*NT);
    for(let dt=-4;dt<=4;dt++){const t=((t0+dt)%NT+NT)%NT,r=Math.round(x*cs[t]+y*sn[t])+diag;acc[t*NR+r]++;}
  }
  const peaks=[];
  for(let t=0;t<NT;t++)for(let r=0;r<NR;r++){
    const v=acc[t*NR+r];if(v<25)continue;let mx=true;
    for(let a=-3;a<=3&&mx;a++)for(let b=-4;b<=4;b++){if(!a&&!b)continue;let tt=t+a,rr=r+b;if(tt<0){tt+=NT;rr=NR-1-rr;}else if(tt>=NT){tt-=NT;rr=NR-1-rr;}
      if(rr<0||rr>=NR)continue;const u=acc[tt*NR+rr];if(u>v||(u===v&&(a<0||(a===0&&b<0)))){mx=false;break;}}
    if(mx)peaks.push({th:t*Math.PI/NT,rho:r-diag,votes:v});
  }
  return peaks.sort((a,b)=>b.votes-a.votes).slice(0,60);
}
function inter(a,b){
  const c1=Math.cos(a.th),s1=Math.sin(a.th),c2=Math.cos(b.th),s2=Math.sin(b.th),d=c1*s2-s1*c2;if(Math.abs(d)<1e-6)return null;
  return[(a.rho*s2-b.rho*s1)/d,(c1*b.rho-c2*a.rho)/d];
}
const sd=(l,x,y)=>x*Math.cos(l.th)+y*Math.sin(l.th)-l.rho;
const angDiff=(a,b)=>{let d=Math.abs(a-b)%Math.PI;return Math.min(d,Math.PI-d);};
/** Fracción de los lados del cuadrilátero q (4 esquinas) donde hay un borde fuerte perpendicular. */
function support(E,w,h,q,thr){
  let ok=0,n=0;
  for(let i=0;i<4;i++){
    const[a,b]=[q[i],q[(i+1)%4]],L=Math.hypot(b[0]-a[0],b[1]-a[1]);if(L<8)return 0;
    const nx=-(b[1]-a[1])/L,ny=(b[0]-a[0])/L,K=Math.max(12,Math.round(L/4));
    for(let k=0;k<K;k++){const t=(k+.5)/K,x=a[0]+(b[0]-a[0])*t,y=a[1]+(b[1]-a[1])*t;
      if(x<2||y<2||x>w-3||y>h-3)continue;n++;
      let best=0;for(let d=-2;d<=2;d++){const X=x+nx*d,Y=y+ny*d,g=Math.abs(samp(E.gx,w,h,X,Y)*nx+samp(E.gy,w,h,X,Y)*ny);if(g>best)best=g;}
      if(best>thr)ok++;}
  }
  return n>=24?ok/n:0;
}
/** Ajusta cada lado a los puntos de borde cercanos (±4 px) y devuelve las esquinas afinadas. */
function refine(E,w,h,q,thr){
  const lines=[];
  for(let i=0;i<4;i++){
    const[a,b]=[q[i],q[(i+1)%4]],L=Math.hypot(b[0]-a[0],b[1]-a[1]),nx=-(b[1]-a[1])/L,ny=(b[0]-a[0])/L,K=Math.max(16,Math.round(L/2)),P=[];
    for(let k=0;k<K;k++){const t=.05+.9*(k+.5)/K,x=a[0]+(b[0]-a[0])*t,y=a[1]+(b[1]-a[1])*t;
      if(x<3||y<3||x>w-4||y>h-4)continue;let bd=0,bv=0;
      for(let d=-4;d<=4;d+=.5){const X=x+nx*d,Y=y+ny*d,g=Math.abs(samp(E.gx,w,h,X,Y)*nx+samp(E.gy,w,h,X,Y)*ny);if(g>bv){bv=g;bd=d;}}
      if(bv>thr)P.push([x+nx*bd,y+ny*bd]);}
    if(P.length<6)return q;
    // recta por mínimos cuadrados (dirección principal)
    let mx=0,my=0;for(const p of P){mx+=p[0];my+=p[1];}mx/=P.length;my/=P.length;
    let sxx=0,sxy=0,syy=0;for(const p of P){const dx=p[0]-mx,dy=p[1]-my;sxx+=dx*dx;sxy+=dx*dy;syy+=dy*dy;}
    const th=.5*Math.atan2(2*sxy,sxx-syy)+Math.PI/2;lines.push({th,rho:mx*Math.cos(th)+my*Math.sin(th)});
  }
  const out=[];for(let i=0;i<4;i++){const p=inter(lines[(i+3)%4],lines[i]);if(!p)return q;out.push(p);}
  for(let i=0;i<4;i++)if(Math.hypot(out[i][0]-q[i][0],out[i][1]-q[i][1])>12)return q;
  return out;
}
function convexContains(q,x,y){
  let s=0;for(let i=0;i<4;i++){const[a,b]=[q[i],q[(i+1)%4]],c=(b[0]-a[0])*(y-a[1])-(b[1]-a[1])*(x-a[0]);if(!s)s=Math.sign(c);else if(Math.sign(c)!==s)return false;}
  return true;
}
function convex(q){let s=0;for(let i=0;i<4;i++){const[a,b,c]=[q[i],q[(i+1)%4],q[(i+2)%4]],cr=(b[0]-a[0])*(c[1]-b[1])-(b[1]-a[1])*(c[0]-b[0]);if(!cr)return false;if(!s)s=Math.sign(cr);else if(Math.sign(cr)!==s)return false;}return true;}
/** Ordena las esquinas: sup-izq, sup-der, inf-der, inf-izq. */
function order(q){
  const cx=q.reduce((s,p)=>s+p[0],0)/4,cy=q.reduce((s,p)=>s+p[1],0)/4;
  const r=q.slice().sort((a,b)=>Math.atan2(a[1]-cy,a[0]-cx)-Math.atan2(b[1]-cy,b[0]-cx));
  let k=0,best=1e18;for(let i=0;i<4;i++){const v=r[i][0]+r[i][1];if(v<best){best=v;k=i;}}
  return[0,1,2,3].map(i=>r[(k+i)%4]);
}

/**
 * Cuadriláteros candidatos alrededor del punto (tx,ty). gray: Uint8Array w×h.
 * Devuelve hasta 8 candidatos ordenados (el mejor primero), cada uno {q:[[x,y]×4], score}.
 */
export function findQuads(gray,w,h,tx,ty){
  const I=new Float32Array(w*h);for(let i=0;i<w*h;i++)I[i]=gray[i];
  // suavizado leve para que el ruido no cuente como borde
  const B=new Float32Array(w*h);
  for(let y=1;y<h-1;y++)for(let x=1;x<w-1;x++){const o=y*w+x;B[o]=(I[o-w-1]+2*I[o-w]+I[o-w+1]+2*I[o-1]+4*I[o]+2*I[o+1]+I[o+w-1]+2*I[o+w]+I[o+w+1])/16;}
  const E=sobel(B,w,h);
  const ms=Array.from(E.mag).filter(v=>v>0).sort((a,b)=>a-b),thr=Math.max(40,ms.length?ms[Math.floor(ms.length*.9)]:40);
  const L=hough(E,w,h,thr).slice(0,28); // las líneas más largas: los lados de un marco o una pantalla
  // para cada orientación, las líneas más cercanas a cada lado del punto
  const pairs=[];
  for(let i=0;i<L.length;i++)for(let j=i+1;j<L.length;j++){
    const a=L[i],b=L[j];if(angDiff(a.th,b.th)>.6)continue;
    const da=sd(a,tx,ty),db=sd(b,tx,ty);
    // a y b a lados opuestos del punto (con orientaciones casi opuestas, el signo se invierte)
    const same=Math.abs(a.th-b.th)<Math.PI/2;
    if(same?Math.sign(da)===Math.sign(db):Math.sign(da)!==Math.sign(db))continue;
    if(Math.abs(da)<4||Math.abs(db)<4)continue;
    pairs.push({a,b,th:a.th,dist:Math.abs(da)+Math.abs(db)});
  }
  const P=pairs,cands=[],area0=w*h;
  for(let i=0;i<P.length;i++)for(let j=i+1;j<P.length;j++){
    const p=P[i],r=P[j];if(angDiff(p.th,r.th)<.5)continue;
    const c=[inter(p.a,r.a),inter(r.a,p.b),inter(p.b,r.b),inter(r.b,p.a)];
    if(c.some(v=>!v||!isFinite(v[0])||!isFinite(v[1])||v[0]<-w*.15||v[1]<-h*.15||v[0]>w*1.15||v[1]>h*1.15))continue;
    const q=order(c);if(!convex(q)||!convexContains(q,tx,ty))continue;
    let area=0;for(let k=0;k<4;k++){const[a,b]=[q[k],q[(k+1)%4]];area+=a[0]*b[1]-b[0]*a[1];}area=Math.abs(area)/2;
    if(area<area0*.005)continue;
    const s=support(E,w,h,q,thr*.6);if(s<.45)continue;
    cands.push({q,score:s*(.75+.25*Math.sqrt(area/area0)),area});
  }
  cands.sort((x,y)=>y.score-x.score);
  // Casi siempre se quiere cubrir lo de adentro (el lienzo de un cuadro, la pantalla de un televisor), no el
  // marco: si dentro del mejor hay otro casi igual de bien marcado y de tamaño parecido, va primero el más chico.
  if(cands.length){
    const top=cands[0],inside=c=>c.q.every(p=>convexContains(top.q,p[0],p[1]));
    const inner=cands.filter(c=>c!==top&&c.score>=top.score*.92&&c.area>=top.area*.45&&c.area<top.area*.97&&inside(c)).sort((a,b)=>a.area-b.area)[0];
    if(inner){cands.splice(cands.indexOf(inner),1);cands.unshift(inner);}
  }
  // sin repetidos (cuadriláteros casi iguales)
  const out=[];
  for(const c of cands){
    if(out.some(o=>o.q.every((p,k)=>Math.hypot(p[0]-c.q[k][0],p[1]-c.q[k][1])<Math.max(6,Math.sqrt(c.area)*.05))))continue;
    out.push({q:refine(E,w,h,c.q,thr*.6),score:c.score});if(out.length>=8)break;
  }
  return out;
}
