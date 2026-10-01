/* Oblicua · seguimiento de una superficie plana a lo largo de un video.
   Sin dependencias y todo en el navegador:
   - Lucas-Kanade piramidal sigue puntos de la superficie de un fotograma al siguiente.
   - Cada fotograma se ajusta además contra la imagen de referencia (el fotograma donde marcaste
     las esquinas), así el letrero no se va desviando con el tiempo.
   - Descriptores binarios estilo ORB vuelven a encontrar la superficie cuando sale del encuadre
     y entra de nuevo.
   - El video se recorre hacia adelante y hacia atrás y se combinan las dos pasadas.
   Coordenadas internas: píxeles del fotograma de análisis, con el centro del píxel i en i. */

const LEVELS=4;      // niveles de la pirámide para el flujo óptico
const R=4;           // radio de la ventana de Lucas-Kanade (9×9)
const ITERS=12;
const MARGIN=.12;    // fracción alrededor de la superficie que también se usa (la pared es el mismo plano)
const MIN_EIG=1.5;   // textura mínima en la ventana para poder seguir un punto
const MAX_PTS=120;

/* ---------- números aleatorios reproducibles ---------- */
function rng(seed){return()=>{seed|=0;seed=seed+0x6D2B79F5|0;let t=Math.imul(seed^seed>>>15,1|seed);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}

/* ---------- homografías: (u,v) del cuadrado unidad → píxel, fila mayor ---------- */
export function hApply(H,u,v){const w=H[6]*u+H[7]*v+H[8];return[(H[0]*u+H[1]*v+H[2])/w,(H[3]*u+H[4]*v+H[5])/w];}
export function hInv(m){
  const[a,b,c,d,e,f,g,h,i]=m,A=e*i-f*h,B=-(d*i-f*g),C=d*h-e*g,det=a*A+b*B+c*C;
  return Float64Array.of(A/det,-(b*i-c*h)/det,(b*f-c*e)/det,B/det,(a*i-c*g)/det,-(a*f-c*d)/det,C/det,-(a*h-b*g)/det,(a*e-b*d)/det);
}
export function hMul(A,B){const o=new Float64Array(9);for(let r=0;r<3;r++)for(let c=0;c<3;c++)o[r*3+c]=A[r*3]*B[c]+A[r*3+1]*B[3+c]+A[r*3+2]*B[6+c];return o;}
/** Esquinas [x0,y0,…,x3,y3] (sup-izq, sup-der, inf-der, inf-izq) → homografía. */
export function hFromQuad(q){
  const[x0,y0,x1,y1,x2,y2,x3,y3]=q,dx1=x1-x2,dx2=x3-x2,dx3=x0-x1+x2-x3,dy1=y1-y2,dy2=y3-y2,dy3=y0-y1+y2-y3;
  let a,b,c,d,e,f,g,h;
  if(Math.abs(dx3)<1e-9&&Math.abs(dy3)<1e-9){a=x1-x0;b=x2-x1;c=x0;d=y1-y0;e=y2-y1;f=y0;g=0;h=0;}
  else{const det=dx1*dy2-dx2*dy1;g=(dx3*dy2-dx2*dy3)/det;h=(dx1*dy3-dx3*dy1)/det;a=x1-x0+g*x1;b=x3-x0+h*x3;c=x0;d=y1-y0+g*y1;e=y3-y0+h*y3;f=y0;}
  return Float64Array.of(a,b,c,d,e,f,g,h,1);
}
export function hQuad(H){return[...hApply(H,0,0),...hApply(H,1,0),...hApply(H,1,1),...hApply(H,0,1)];}

function solve8(A,b){
  const M=new Float64Array(72);
  for(let r=0;r<8;r++){for(let c=0;c<8;c++)M[r*9+c]=A[r*8+c];M[r*9+8]=b[r];}
  for(let c=0;c<8;c++){
    let p=c,mx=Math.abs(M[c*9+c]);
    for(let r=c+1;r<8;r++){const v=Math.abs(M[r*9+c]);if(v>mx){mx=v;p=r;}}
    if(mx<1e-12)return null;
    if(p!==c)for(let k=c;k<9;k++){const t=M[c*9+k];M[c*9+k]=M[p*9+k];M[p*9+k]=t;}
    for(let r=c+1;r<8;r++){const f=M[r*9+c]/M[c*9+c];if(f)for(let k=c;k<9;k++)M[r*9+k]-=f*M[c*9+k];}
  }
  const x=new Float64Array(8);
  for(let r=7;r>=0;r--){let s=M[r*9+8];for(let k=r+1;k<8;k++)s-=M[r*9+k]*x[k];x[r]=s/M[r*9+r];}
  return x;
}
/** Homografía por mínimos cuadrados (DLT normalizado) de S[idx] → D[idx]. */
function fitH(S,D,idx,n){
  let mx=0,my=0,nx=0,ny=0;
  for(let k=0;k<n;k++){const i=idx[k];mx+=S[2*i];my+=S[2*i+1];nx+=D[2*i];ny+=D[2*i+1];}
  mx/=n;my/=n;nx/=n;ny/=n;
  let ds=0,dd=0;
  for(let k=0;k<n;k++){const i=idx[k];ds+=Math.hypot(S[2*i]-mx,S[2*i+1]-my);dd+=Math.hypot(D[2*i]-nx,D[2*i+1]-ny);}
  if(ds<1e-12||dd<1e-12)return null;
  const ss=Math.SQRT2*n/ds,sd=Math.SQRT2*n/dd,A=new Float64Array(64),b=new Float64Array(8),r1=new Float64Array(8),r2=new Float64Array(8);
  for(let k=0;k<n;k++){
    const i=idx[k],u=(S[2*i]-mx)*ss,v=(S[2*i+1]-my)*ss,x=(D[2*i]-nx)*sd,y=(D[2*i+1]-ny)*sd;
    r1[0]=u;r1[1]=v;r1[2]=1;r1[3]=0;r1[4]=0;r1[5]=0;r1[6]=-u*x;r1[7]=-v*x;
    r2[0]=0;r2[1]=0;r2[2]=0;r2[3]=u;r2[4]=v;r2[5]=1;r2[6]=-u*y;r2[7]=-v*y;
    for(let p=0;p<8;p++){const a1=r1[p],a2=r2[p];b[p]+=a1*x+a2*y;for(let q=p;q<8;q++)A[p*8+q]+=a1*r1[q]+a2*r2[q];}
  }
  for(let p=0;p<8;p++)for(let q=0;q<p;q++)A[p*8+q]=A[q*8+p];
  const h=solve8(A,b);if(!h)return null;
  const Hn=Float64Array.of(h[0],h[1],h[2],h[3],h[4],h[5],h[6],h[7],1);
  const Ts=Float64Array.of(ss,0,-ss*mx,0,ss,-ss*my,0,0,1),Td=Float64Array.of(1/sd,0,nx,0,1/sd,ny,0,0,1);
  const H=hMul(Td,hMul(Hn,Ts));
  return H.every(Number.isFinite)?H:null;
}
function countIn(H,S,D,n,th2,mark){
  let c=0;
  for(let i=0;i<n;i++){
    const u=S[2*i],v=S[2*i+1],w=H[6]*u+H[7]*v+H[8];
    const ex=(H[0]*u+H[1]*v+H[2])/w-D[2*i],ey=(H[3]*u+H[4]*v+H[5])/w-D[2*i+1],ok=ex*ex+ey*ey<th2;
    if(mark)mark[i]=ok?1:0; if(ok)c++;
  }
  return c;
}
/** RANSAC: S (u,v) → D (píxeles). Devuelve {H,cnt,inl} o null. */
function ransac(S,D,n,th,iters,minInl,rand){
  if(n<Math.max(4,minInl))return null;
  const th2=th*th,smp=new Int32Array(4);let best=0,bestH=null,N=iters;
  for(let it=0;it<N;it++){
    for(let k=0;k<4;k++){let r;do{r=rand()*n|0;}while(smp.slice(0,k).includes(r));smp[k]=r;}
    const H=fitH(S,D,smp,4);if(!H)continue;
    const c=countIn(H,S,D,n,th2,null);
    if(c>best){best=c;bestH=H;const e=c/n;if(e>=1)break;const need=Math.log(1-.995)/Math.log(1-e**4);if(need+10<N)N=Math.max(it+1,need+10|0);}
  }
  if(best<minInl)return null;
  const inl=new Uint8Array(n);let H=bestH;
  for(let r=0;r<2;r++){
    countIn(H,S,D,n,th2,inl);
    const idx=[];for(let i=0;i<n;i++)if(inl[i])idx.push(i);
    if(idx.length<4)break;
    const H2=fitH(S,D,idx,idx.length);if(!H2)break;H=H2;
  }
  const cnt=countIn(H,S,D,n,th2,inl);
  return cnt>=minInl?{H,cnt,inl}:null;
}

/* ---------- imágenes ---------- */
function samp(I,w,h,x,y){
  if(x<0)x=0;else if(x>w-1.001)x=w-1.001;
  if(y<0)y=0;else if(y>h-1.001)y=h-1.001;
  const ix=x|0,iy=y|0,fx=x-ix,fy=y-iy,o=iy*w+ix,a=I[o]+(I[o+1]-I[o])*fx,b=I[o+w]+(I[o+w+1]-I[o+w])*fx;
  return a+(b-a)*fy;
}
function grads(L){
  const{w,h,I}=L,X=new Float32Array(w*h),Y=new Float32Array(w*h);
  for(let y=1;y<h-1;y++)for(let x=1;x<w-1;x++){
    const o=y*w+x,a=I[o-w-1],b=I[o-w],c=I[o-w+1],d=I[o-1],f=I[o+1],g=I[o+w-1],k=I[o+w],l=I[o+w+1];
    X[o]=(3*(c-a)+10*(f-d)+3*(l-g))/32; Y[o]=(3*(g-a)+10*(k-b)+3*(l-c))/32;
  }
  L.Ix=X;L.Iy=Y;return L;
}
/** Desenfoque [1 2 1] y submuestreo: el píxel i del nivel siguiente cae sobre el píxel 2i. */
function down(L){
  const{w,h,I}=L,w2=(w+1)>>1,h2=(h+1)>>1,O=new Float32Array(w2*h2);
  for(let y=0;y<h2;y++){const yc=2*y,ym=yc>0?yc-1:0,yp=yc+1<h?yc+1:yc;
    const a=ym*w,b=yc*w,c=yp*w;
    for(let x=0;x<w2;x++){const xc=2*x,xm=xc>0?xc-1:0,xp=xc+1<w?xc+1:xc;
      O[y*w2+x]=(I[a+xm]+2*I[a+xc]+I[a+xp]+2*(I[b+xm]+2*I[b+xc]+I[b+xp])+I[c+xm]+2*I[c+xc]+I[c+xp])/16;}}
  return{w:w2,h:h2,I:O};
}
function pyramid(I,w,h,n){const P=[grads({w,h,I})];while(P.length<n&&P[P.length-1].w>24&&P[P.length-1].h>24)P.push(grads(down(P[P.length-1])));return P;}
function blur(I,w,h,k){ // separable, k = pesos simétricos (centro primero)
  const T=new Float32Array(w*h),O=new Float32Array(w*h),r=k.length-1;let s=k[0];for(let i=1;i<=r;i++)s+=2*k[i];
  for(let y=0;y<h;y++)for(let x=0;x<w;x++){let a=I[y*w+x]*k[0];for(let i=1;i<=r;i++)a+=k[i]*(I[y*w+(x-i<0?0:x-i)]+I[y*w+(x+i>=w?w-1:x+i)]);T[y*w+x]=a/s;}
  for(let y=0;y<h;y++)for(let x=0;x<w;x++){let a=T[y*w+x]*k[0];for(let i=1;i<=r;i++)a+=k[i]*(T[(y-i<0?0:y-i)*w+x]+T[(y+i>=h?h-1:y+i)*w+x]);O[y*w+x]=a/s;}
  return O;
}
function resize(I,w,h,s){
  const w2=Math.floor((w-1)*s)+1,h2=Math.floor((h-1)*s)+1,O=new Float32Array(w2*h2);
  for(let y=0;y<h2;y++)for(let x=0;x<w2;x++)O[y*w2+x]=samp(I,w,h,x/s,y/s);
  return{w:w2,h:h2,I:O};
}

/* ---------- Lucas-Kanade piramidal (con compensación de brillo) ---------- */
const WN=(2*R+1)**2,wI=new Float32Array(WN),wX=new Float32Array(WN),wY=new Float32Array(WN),wJ=new Float32Array(WN);
/** Sigue n puntos de la pirámide A (posiciones PA) a la pirámide B (estimación inicial PB). */
function lk(A,B,PA,PB,n,out,ok){
  const top=Math.min(A.length,B.length)-1,B0=B[0];
  for(let i=0;i<n;i++){
    let bx=PB[2*i]/(1<<top),by=PB[2*i+1]/(1<<top),good=true,err=0;
    for(let L=top;L>=0;L--){
      const a=A[L],b=B[L],s=1/(1<<L),ax=PA[2*i]*s,ay=PA[2*i+1]*s;
      let gxx=0,gxy=0,gyy=0,mI=0,k=0;
      for(let dy=-R;dy<=R;dy++)for(let dx=-R;dx<=R;dx++,k++){
        const x=ax+dx,y=ay+dy,vI=samp(a.I,a.w,a.h,x,y),vX=samp(a.Ix,a.w,a.h,x,y),vY=samp(a.Iy,a.w,a.h,x,y);
        wI[k]=vI;wX[k]=vX;wY[k]=vY;mI+=vI;gxx+=vX*vX;gxy+=vX*vY;gyy+=vY*vY;
      }
      mI/=WN;
      const det=gxx*gyy-gxy*gxy,eig=(gxx+gyy-Math.sqrt((gxx-gyy)**2+4*gxy*gxy))/2/WN;
      if(eig<MIN_EIG||det<=0){ if(L===0){good=false;break;} bx*=2;by*=2;continue; }
      for(let it=0;it<ITERS;it++){
        let mJ=0;k=0;
        for(let dy=-R;dy<=R;dy++)for(let dx=-R;dx<=R;dx++,k++){const v=samp(b.I,b.w,b.h,bx+dx,by+dy);wJ[k]=v;mJ+=v;}
        mJ/=WN;let sx=0,sy=0;err=0;
        for(k=0;k<WN;k++){const d=(wI[k]-mI)-(wJ[k]-mJ);sx+=d*wX[k];sy+=d*wY[k];err+=Math.abs(d);}
        const nx=(gyy*sx-gxy*sy)/det,ny=(gxx*sy-gxy*sx)/det;
        bx+=nx;by+=ny;
        if(nx*nx+ny*ny<1e-4)break;
      }
      if(L>0){bx*=2;by*=2;}
    }
    out[2*i]=bx;out[2*i+1]=by;
    ok[i]=good&&err/WN<28&&bx>=R&&by>=R&&bx<=B0.w-1-R&&by<=B0.h-1-R?1:0;
  }
}

/* ---------- esquinas (Shi-Tomasi) ---------- */
/** Mejores esquinas de un nivel dentro de la caja reg (y de reg.test si existe). */
function corners(L,reg,maxN,cell,avoid){
  const{w,Ix,Iy}=L,x0=Math.max(reg.x0,3),y0=Math.max(reg.y0,3),x1=Math.min(reg.x1,L.w-4),y1=Math.min(reg.y1,L.h-4);
  if(x1-x0<4||y1-y0<4)return[];
  const bw=x1-x0+5,bh=y1-y0+5,xx=new Float32Array(bw*bh),xy=new Float32Array(bw*bh),yy=new Float32Array(bw*bh);
  for(let y=0;y<bh;y++)for(let x=0;x<bw;x++){const o=(y+y0-2)*w+x+x0-2,gx=Ix[o],gy=Iy[o],j=y*bw+x;xx[j]=gx*gx;xy[j]=gx*gy;yy[j]=gy*gy;}
  const box=a=>{const t=new Float32Array(bw*bh);for(let y=0;y<bh;y++)for(let x=2;x<bw-2;x++){const j=y*bw+x;t[j]=a[j-2]+a[j-1]+a[j]+a[j+1]+a[j+2];}
    const o=new Float32Array(bw*bh);for(let y=2;y<bh-2;y++)for(let x=0;x<bw;x++){const j=y*bw+x;o[j]=t[j-2*bw]+t[j-bw]+t[j]+t[j+bw]+t[j+2*bw];}return o;};
  const A=box(xx),Bm=box(xy),C=box(yy),E=new Float32Array(bw*bh);let mx=0;
  for(let y=2;y<bh-2;y++)for(let x=2;x<bw-2;x++){const j=y*bw+x,a=A[j],b=Bm[j],c=C[j],e=(a+c-Math.sqrt((a-c)**2+4*b*b))/50;E[j]=e;if(e>mx)mx=e;}
  const thr=Math.max(MIN_EIG*2,mx*.02),cells=new Map();
  for(let y=3;y<bh-3;y++)for(let x=3;x<bw-3;x++){
    const j=y*bw+x,e=E[j];if(e<thr)continue;
    if(e<E[j-1]||e<=E[j+1]||e<E[j-bw]||e<=E[j+bw]||e<E[j-bw-1]||e<E[j-bw+1]||e<=E[j+bw-1]||e<=E[j+bw+1])continue;
    const px=x+x0-2,py=y+y0-2;if(reg.test&&!reg.test(px,py))continue;
    const ck=((py/cell)|0)*4096+((px/cell)|0);if(avoid&&avoid.has(ck))continue;
    const c=cells.get(ck);if(!c||c[2]<e)cells.set(ck,[px,py,e]);
  }
  return[...cells.values()].sort((a,b)=>b[2]-a[2]).slice(0,maxN);
}

/** Zona de la superficie (más el margen) en el fotograma, recortada a la imagen. */
function region(H,w,h,border,m=MARGIN){
  const Hi=hInv(H),c=hApply(H,.5,.5),sg=Math.sign(Hi[6]*c[0]+Hi[7]*c[1]+Hi[8]);
  let X0=1e9,Y0=1e9,X1=-1e9,Y1=-1e9;
  for(const[u,v]of[[-m,-m],[1+m,-m],[1+m,1+m],[-m,1+m]]){const p=hApply(H,u,v);if(!isFinite(p[0])||!isFinite(p[1]))return null;X0=Math.min(X0,p[0]);Y0=Math.min(Y0,p[1]);X1=Math.max(X1,p[0]);Y1=Math.max(Y1,p[1]);}
  const x0=Math.max(border,Math.floor(X0)),y0=Math.max(border,Math.floor(Y0)),x1=Math.min(w-1-border,Math.ceil(X1)),y1=Math.min(h-1-border,Math.ceil(Y1));
  if(x1-x0<6||y1-y0<6)return null;
  const test=(x,y)=>{const d=Hi[6]*x+Hi[7]*y+Hi[8];if(d*sg<=0)return false;const u=(Hi[0]*x+Hi[1]*y+Hi[2])/d,v=(Hi[3]*x+Hi[4]*y+Hi[5])/d;return u>=-m&&u<=1+m&&v>=-m&&v<=1+m;};
  return{x0,y0,x1,y1,test,Hi};
}

/** ¿La homografía describe una superficie creíble en un fotograma w×h? ori: sentido de giro de las esquinas. */
function sane(H,w,h,ori){
  if(!H)return false;
  const ds=[[0,0],[1,0],[1,1],[0,1]].map(([u,v])=>H[6]*u+H[7]*v+H[8]),s0=Math.sign(H[6]*.5+H[7]*.5+H[8]);
  if(ds.some(d=>d*s0<=0))return false;
  const ad=ds.map(Math.abs);if(Math.max(...ad)/Math.min(...ad)>12)return false;
  const q=hQuad(H);if(!q.every(Number.isFinite))return false;
  let area=0;
  for(let i=0;i<4;i++){
    const a=i*2,b=(i+1)%4*2,c=(i+2)%4*2,cr=(q[b]-q[a])*(q[c+1]-q[b+1])-(q[b+1]-q[a+1])*(q[c]-q[b]);
    if(Math.sign(cr)!==ori)return false;
    area+=q[a]*q[b+1]-q[b]*q[a+1];
    if(Math.hypot(q[b]-q[a],q[b+1]-q[a+1])<3)return false;
  }
  area=Math.abs(area)/2;if(area<w*h*.0004||area>w*h*40)return false;
  const xs=[q[0],q[2],q[4],q[6]],ys=[q[1],q[3],q[5],q[7]];
  return Math.max(...xs)>0&&Math.min(...xs)<w&&Math.max(...ys)>0&&Math.min(...ys)<h;
}
function orient(H){const q=hQuad(H);return Math.sign((q[2]-q[0])*(q[5]-q[3])-(q[3]-q[1])*(q[4]-q[2]));}
function jump(A,B,w,h){const a=hQuad(A),b=hQuad(B);let s=0;for(let i=0;i<4;i++)s+=Math.hypot(a[2*i]-b[2*i],a[2*i+1]-b[2*i+1]);return s/4/Math.hypot(w,h);}

/* ---------- descriptores binarios (estilo ORB) ---------- */
const NB=30,PAT=(()=>{
  const r=rng(1234),g=()=>Math.sqrt(-2*Math.log(r()+1e-12))*Math.cos(2*Math.PI*r()),P=new Float32Array(1024);
  for(let i=0;i<1024;i++){let v;do v=g()*5.2;while(Math.abs(v)>11.5);P[i]=v;}
  const O=new Int8Array(NB*1024);
  for(let b=0;b<NB;b++){const a=b*2*Math.PI/NB,c=Math.cos(a),s=Math.sin(a);
    for(let i=0;i<1024;i+=2){O[b*1024+i]=Math.round(c*P[i]-s*P[i+1]);O[b*1024+i+1]=Math.round(s*P[i]+c*P[i+1]);}}
  return O;
})();
const CIRC=(()=>{const o=[];for(let y=-12;y<=12;y++)for(let x=-12;x<=12;x++)if(x*x+y*y<=144)o.push(x,y);return Int8Array.from(o);})();
const FB=16; // borde mínimo para un descriptor
function popc(v){v-=(v>>>1)&0x55555555;v=(v&0x33333333)+((v>>>2)&0x33333333);return(((v+(v>>>4))&0x0F0F0F0F)*0x01010101)>>>24;}
/** Puntos con descriptor en 4 escalas. reg (opcional, en coordenadas del nivel 0) limita la zona. */
function features(pyr,reg,maxN){
  const L0=pyr[0],L1=pyr[1]||null,sc=[];
  sc.push({s:1,L:L0});
  sc.push({s:Math.SQRT1_2,L:grads(resize(blur(L0.I,L0.w,L0.h,[2,1]),L0.w,L0.h,Math.SQRT1_2))});
  if(L1){sc.push({s:.5,L:L1});sc.push({s:.5*Math.SQRT1_2,L:grads(resize(blur(L1.I,L1.w,L1.h,[2,1]),L1.w,L1.h,Math.SQRT1_2))});}
  const tot=sc.reduce((a,c)=>a+c.s*c.s,0),X=[],Y=[],Dd=[];
  for(const{s,L}of sc){
    if(L.w<2*FB+8||L.h<2*FB+8)continue;
    let r={x0:FB,y0:FB,x1:L.w-1-FB,y1:L.h-1-FB};
    if(reg){r={x0:Math.max(FB,Math.floor(reg.x0*s)),y0:Math.max(FB,Math.floor(reg.y0*s)),x1:Math.min(L.w-1-FB,Math.ceil(reg.x1*s)),y1:Math.min(L.h-1-FB,Math.ceil(reg.y1*s)),test:(x,y)=>reg.test(x/s,y/s)};if(r.x1-r.x0<6||r.y1-r.y0<6)continue;}
    const n=Math.max(12,Math.round(maxN*s*s/tot)),cs=corners(L,r,n,Math.max(5,Math.round(Math.sqrt((r.x1-r.x0)*(r.y1-r.y0)/n)*.8)),null);
    if(!cs.length)continue;
    const D=blur(L.I,L.w,L.h,[6,4,1]),w=L.w;
    for(const[x,y]of cs){
      let m01=0,m10=0;const o=y*w+x;
      for(let i=0;i<CIRC.length;i+=2){const v=D[o+CIRC[i+1]*w+CIRC[i]];m10+=CIRC[i]*v;m01+=CIRC[i+1]*v;}
      let b=Math.round(Math.atan2(m01,m10)/(2*Math.PI)*NB);b=((b%NB)+NB)%NB;
      const P=b*1024;
      for(let k=0;k<8;k++){let word=0;
        for(let j=0;j<32;j++){const q=P+(k*32+j)*4;if(D[o+PAT[q+1]*w+PAT[q]]<D[o+PAT[q+3]*w+PAT[q+2]])word|=1<<j;}
        Dd.push(word>>>0);}
      X.push(x/s);Y.push(y/s);
    }
  }
  return{x:Float32Array.from(X),y:Float32Array.from(Y),d:Uint32Array.from(Dd),n:X.length};
}
function matchF(F,G){ // F: fotograma, G: referencia → pares [i_ref, i_fotograma]
  const out=[],a=F.d,b=G.d;
  for(let i=0;i<F.n;i++){
    let b1=999,b2=999,bj=-1;const o=i*8;
    for(let j=0;j<G.n;j++){const p=j*8;
      const d=popc(a[o]^b[p])+popc(a[o+1]^b[p+1])+popc(a[o+2]^b[p+2])+popc(a[o+3]^b[p+3])+popc(a[o+4]^b[p+4])+popc(a[o+5]^b[p+5])+popc(a[o+6]^b[p+6])+popc(a[o+7]^b[p+7]);
      if(d<b1){b2=b1;b1=d;bj=j;}else if(d<b2)b2=d;}
    if(bj>=0&&b1<=72&&b1<.82*b2)out.push(bj,i);
  }
  return out;
}

/* ---------- referencia: un fotograma donde se marcaron las esquinas ---------- */
/** gray: Uint8Array w×h. quad: esquinas en píxeles de análisis. */
export function buildRef(gray,w,h,quad){
  const I=new Float32Array(w*h);for(let i=0;i<w*h;i++)I[i]=gray[i];
  const pyr=pyramid(I,w,h,LEVELS),H=hFromQuad(quad),ori=orient(H);
  const reg=region(H,w,h,FB);
  const F=reg?features(pyr,reg,450):{x:new Float32Array(0),y:new Float32Array(0),d:new Uint32Array(0),n:0};
  const Hi=hInv(H),u=new Float32Array(F.n),v=new Float32Array(F.n);
  for(let i=0;i<F.n;i++){const p=hApply(Hi,F.x[i],F.y[i]);u[i]=p[0];v[i]=p[1];}
  const reg2=region(H,w,h,R+2),pts=[];
  if(reg2)for(const[x,y]of corners(pyr[0],reg2,260,7,null)){const p=hApply(Hi,x,y);pts.push(p[0],p[1]);}
  return{pyr,H,ori,F:{...F,u,v},pts:Float32Array.from(pts),w,h};
}

/** Afina H contra la imagen de referencia: así no se acumula error de un fotograma a otro. */
function anchor(pyr,H,ref,rand){
  const L0=pyr[0],w=L0.w,h=L0.h,reg=region(H,w,h,R+2);if(!reg)return null;
  const{x0,y0,x1,y1}=reg,bw=x1-x0+1,bh=y1-y0+1;if(bw<20||bh<20)return null;
  const M=hMul(ref.H,reg.Hi),c=hApply(H,.5,.5),p0=hApply(M,c[0],c[1]),p1=hApply(M,c[0]+1,c[1]),p2=hApply(M,c[0],c[1]+1);
  const sc=Math.sqrt(Math.abs((p1[0]-p0[0])*(p2[1]-p0[1])-(p1[1]-p0[1])*(p2[0]-p0[0])));
  let rl=0;while(rl<ref.pyr.length-1&&sc/(1<<rl)>1.6)rl++;
  const RL=ref.pyr[rl],k=1/(1<<rl),I=new Float32Array(bw*bh);
  for(let y=0;y<bh;y++){const Y=y+y0;for(let x=0;x<bw;x++){const X=x+x0,d=M[6]*X+M[7]*Y+M[8];
    I[y*bw+x]=samp(RL.I,RL.w,RL.h,(M[0]*X+M[1]*Y+M[2])/d*k,(M[3]*X+M[4]*Y+M[5])/d*k);}}
  const A=pyramid(I,bw,bh,3),PA=[],PB=[],U=[];
  for(let i=0;i<ref.pts.length&&U.length<2*MAX_PTS;i+=2){
    const p=hApply(H,ref.pts[i],ref.pts[i+1]);
    if(p[0]<x0+R+1||p[1]<y0+R+1||p[0]>x1-R-1||p[1]>y1-R-1)continue;
    PA.push(p[0]-x0,p[1]-y0);PB.push(p[0],p[1]);U.push(ref.pts[i],ref.pts[i+1]);
  }
  const n=U.length/2;if(n<8)return null;
  const out=new Float32Array(2*n),ok=new Uint8Array(n);
  lk(A,pyr,Float32Array.from(PA),Float32Array.from(PB),n,out,ok);
  const S=[],D=[];for(let i=0;i<n;i++)if(ok[i]){S.push(U[2*i],U[2*i+1]);D.push(out[2*i],out[2*i+1]);}
  const m=S.length/2,res=ransac(Float64Array.from(S),Float64Array.from(D),m,1.5,200,Math.max(8,Math.round(n*.3)),rand);
  return res&&sane(res.H,w,h,ref.ori)?{H:res.H,cnt:res.cnt}:null;
}

/** Busca la superficie en todo el fotograma comparando descriptores con las referencias. */
function detect(pyr,refs,rand){
  const L0=pyr[0],F=features(pyr,null,520);if(F.n<12)return null;
  let best=null;
  refs.forEach((ref,ri)=>{
    if(ref.F.n<12)return;
    const mt=matchF(F,ref.F),m=mt.length/2;if(m<12)return;
    const S=new Float64Array(2*m),D=new Float64Array(2*m);
    for(let k=0;k<m;k++){const j=mt[2*k],i=mt[2*k+1];S[2*k]=ref.F.u[j];S[2*k+1]=ref.F.v[j];D[2*k]=F.x[i];D[2*k+1]=F.y[i];}
    const res=ransac(S,D,m,3.5,1500,Math.max(12,Math.round(m*.12)),rand);
    if(res&&sane(res.H,L0.w,L0.h,ref.ori)&&(!best||res.cnt>best.cnt))best={H:res.H,cnt:res.cnt,ref:ri};
  });
  return best;
}

/** Homografía que lleva la parametrización de la referencia A a la imagen de la referencia B (o null). */
function relate(A,B,rand){
  const mt=matchF(B.F,A.F),m=mt.length/2;if(m<15)return null;
  const S=new Float64Array(2*m),D=new Float64Array(2*m);
  for(let k=0;k<m;k++){const j=mt[2*k],i=mt[2*k+1];S[2*k]=A.F.u[j];S[2*k+1]=A.F.v[j];D[2*k]=B.F.x[i];D[2*k+1]=B.F.y[i];}
  const res=ransac(S,D,m,3,2000,Math.max(15,Math.round(m*.15)),rand);
  if(!res||!sane(res.H,B.w,B.h,A.ori))return null;
  const a=anchor(B.pyr,res.H,A,rand);return a?a.H:res.H;
}

/* ---------- seguidor fotograma a fotograma ---------- */
export class Tracker{
  constructor(w,h,refs){this.w=w;this.h=h;this.refs=refs;this.rand=rng(99);this.H=null;this.Hp=null;this.pts=null;this.prev=null;this.lost=0;this.age=0;this.ref=0;this.coast=[];}
  /** gray: Uint8Array del fotograma. key: {H,ref} si en este fotograma hay esquinas marcadas a mano. */
  step(gray,key){
    const{w,h}=this,I=new Float32Array(w*h);for(let i=0;i<w*h;i++)I[i]=gray[i];
    const pyr=pyramid(I,w,h,LEVELS);let H=null,reset=false,score=0;
    if(key){H=key.H;this.ref=key.ref;reset=true;score=1e3;}
    else if(this.H&&this.pts&&this.prev&&this.pts.length>=16){
      let Hg=this.H;
      if(this.Hp){const P=hMul(hMul(this.H,hInv(this.Hp)),this.H);if(sane(P,w,h,this.refs[this.ref].ori)&&jump(P,this.H,w,h)<.2)Hg=P;}
      const pts=this.pts,n=pts.length/4,PA=new Float32Array(2*n),PB=new Float32Array(2*n),out=new Float32Array(2*n),ok=new Uint8Array(n);
      for(let i=0;i<n;i++){PA[2*i]=pts[4*i];PA[2*i+1]=pts[4*i+1];const p=hApply(Hg,pts[4*i+2],pts[4*i+3]);PB[2*i]=p[0];PB[2*i+1]=p[1];}
      lk(this.prev,pyr,PA,PB,n,out,ok);
      const S=[],D=[];for(let i=0;i<n;i++)if(ok[i]){S.push(pts[4*i+2],pts[4*i+3]);D.push(out[2*i],out[2*i+1]);}
      const m=S.length/2,res=ransac(Float64Array.from(S),Float64Array.from(D),m,2,250,Math.max(8,Math.round(m*.35)),this.rand);
      const ori=this.refs[this.ref].ori;
      if(res&&sane(res.H,w,h,ori)&&jump(res.H,this.H,w,h)<.3){
        const a=anchor(pyr,res.H,this.refs[this.ref],this.rand);
        if(a&&jump(a.H,res.H,w,h)<.04){H=a.H;score=a.cnt;}else{H=res.H;score=res.cnt*.5;}
      }
    }
    // Un instante sin enganche (desenfoque por movimiento, algo que tapa): sigue con la misma velocidad
    // unos fotogramas. Solo se conservan si luego se vuelve a enganchar.
    let coasting=false;
    if(!H&&!key&&this.H&&this.Hp&&this.coast.length<4){
      const P=hMul(hMul(this.H,hInv(this.Hp)),this.H);
      if(sane(P,w,h,this.refs[this.ref].ori)&&jump(P,this.H,w,h)<.2){H=P;coasting=true;}
    }
    if(!H&&!key){
      this.dropCoast();this.lost++;
      if(this.lost<=4||this.lost%3===0){
        const d=detect(pyr,this.refs,this.rand);
        if(d){const a=anchor(pyr,d.H,this.refs[d.ref],this.rand);H=a?a.H:d.H;score=a?a.cnt:d.cnt*.5;this.ref=d.ref;reset=true;}
      }
    }
    const res={H,age:0,score};
    if(H){
      this.lost=0;this.age=reset?0:this.age+1;res.age=this.age;
      if(coasting){this.coast.push(res);this.Hp=this.H;this.H=H;this.pts=this.points(pyr,H,this.pts,true);}
      else{
        this.coast=[];
        this.Hp=reset?null:this.H;this.H=H;this.pts=this.points(pyr,H,reset?null:this.pts);
      }
    }else{this.H=this.Hp=this.pts=null;}
    this.prev=pyr;
    return res;
  }
  dropCoast(){for(const r of this.coast)r.H=null;this.coast=[];}
  /** Puntos a seguir: los que siguen bien (reproyectados con H) más esquinas nuevas dentro de la superficie. */
  points(pyr,H,old,only){
    const{w,h}=this,reg=region(H,w,h,R+2);if(!reg)return null;
    const keep=[],cell=8,avoid=new Set();
    if(old)for(let i=0;i<old.length;i+=4){
      const u=old[i+2],v=old[i+3],p=hApply(H,u,v);
      if(p[0]<reg.x0||p[1]<reg.y0||p[0]>reg.x1||p[1]>reg.y1)continue;
      const ck=((p[1]/cell)|0)*4096+((p[0]/cell)|0);if(avoid.has(ck))continue;
      avoid.add(ck);keep.push(p[0],p[1],u,v);
    }
    if(!only&&keep.length/4<MAX_PTS*.6){
      for(const[x,y]of corners(pyr[0],reg,MAX_PTS-keep.length/4,cell,avoid)){const p=hApply(reg.Hi,x,y);keep.push(x,y,p[0],p[1]);}
    }
    return Float32Array.from(keep);
  }
}

/* ---------- análisis completo de un video ---------- */
const tk=t=>Math.round(t*1e5);
function grayOf(canvas,w,h){
  const d=canvas.getContext('2d').getImageData(0,0,w,h).data,g=new Uint8Array(w*h);
  for(let i=0,j=0;i<g.length;i++,j+=4)g[i]=(d[j]*77+d[j+1]*150+d[j+2]*29)>>8;
  return g;
}
/** Combina las pasadas hacia adelante y hacia atrás, rellena huecos cortos y suaviza. */
export function merge(T,proc,fw,bw,keyT,w,h,keyU=[]){
  const n=T.length,Q=new Float32Array(n*8).fill(NaN),diag=Math.hypot(w,h),fixed=new Uint8Array(n),done=new Uint8Array(n);
  for(let i=0;i<n;i++){
    const k=tk(T[i]);if(!proc.has(k))continue;done[i]=1;
    if(keyT.has(k))fixed[i]=1;
    const f=fw.get(k),b=bw.get(k),qf=f&&f.H?hQuad(f.H):null,qb=b&&b.H?hQuad(b.H):null;
    let q=null;
    if(qf&&qb){
      let dsum=0;for(let c=0;c<4;c++)dsum+=Math.hypot(qf[2*c]-qb[2*c],qf[2*c+1]-qb[2*c+1]);
      if(dsum/4>.04*diag)q=f.score>=b.score?qf:qb;
      else{const wf=1/(1+f.age),wb=1/(1+b.age);q=qf.map((v,j)=>(v*wf+qb[j]*wb)/(wf+wb));}
    }else q=qf||qb;
    if(q)for(let j=0;j<8;j++)Q[i*8+j]=q[j];
  }
  const valid=i=>!Number.isNaN(Q[i*8]);
  // huecos: fotogramas no analizados o perdidos un instante entre dos válidos
  for(let i=0;i<n;i++){
    if(valid(i))continue;
    let a=i-1;while(a>=0&&!valid(a)&&!done[a])a--;
    let b=i+1;while(b<n&&!valid(b)&&!done[b])b++;
    // en los extremos del video (o junto a un hueco) un fotograma sin analizar toma el valor del vecino
    if((a<0||!valid(a))&&b<n&&valid(b)&&!done[i]&&T[b]-T[i]<.1){for(let j=0;j<8;j++)Q[i*8+j]=Q[b*8+j];continue;}
    if((b>=n||!valid(b))&&a>=0&&valid(a)&&!done[i]&&T[i]-T[a]<.1){for(let j=0;j<8;j++)Q[i*8+j]=Q[a*8+j];continue;}
    if(a<0||b>=n||!valid(a)||!valid(b))continue;
    let lostRun=0;for(let j=a+1;j<b;j++)if(done[j])lostRun++;
    if(lostRun>2||T[b]-T[a]>.2)continue;
    const f=(T[i]-T[a])/(T[b]-T[a]);for(let j=0;j<8;j++)Q[i*8+j]=Q[a*8+j]+(Q[b*8+j]-Q[a*8+j])*f;
  }
  // suavizado ligero (Savitzky-Golay: no se atrasa con movimientos rápidos), sin cruzar huecos,
  // sin mover las esquinas marcadas a mano y como mucho 1 px: solo quita el temblor
  const K=[-3,12,17,12,-3],S=Q.slice();
  for(let i=2;i<n-2;i++){
    if(fixed[i])continue;
    let all=true;for(let o=-2;o<=2;o++)if(!valid(i+o)){all=false;break;}
    if(!all)continue;
    for(let j=0;j<8;j++){let s=0;for(let o=-2;o<=2;o++)s+=Q[(i+o)*8+j]*K[o+2];const d=s/35-Q[i*8+j];S[i*8+j]=Q[i*8+j]+Math.max(-1,Math.min(1,d));}
  }
  // Correcciones que mueven las esquinas sobre la superficie (no solo arreglan el seguimiento):
  // cada momento marcado tiene sus esquinas en coordenadas de la superficie; entre dos momentos se interpolan.
  const U0=[0,0,1,0,1,1,0,1],ks=keyU.slice().sort((a,b)=>a.t-b.t);
  if(ks.some(k=>k.U.some((v,j)=>Math.abs(v-U0[j])>1e-4)))for(let i=0;i<n;i++){
    if(!valid(i))continue;
    const t=T[i];let U;
    if(t<=ks[0].t)U=ks[0].U;else if(t>=ks[ks.length-1].t)U=ks[ks.length-1].U;
    else{let a=0;while(ks[a+1].t<t)a++;const A=ks[a],B=ks[a+1],f=(t-A.t)/(B.t-A.t);U=A.U.map((v,j)=>v+(B.U[j]-v)*f);}
    const Hq=hFromQuad(Array.from(S.subarray(i*8,i*8+8)));
    for(let c=0;c<4;c++){const p=hApply(Hq,U[2*c],U[2*c+1]);S[i*8+2*c]=p[0];S[i*8+2*c+1]=p[1];}
  }
  // a coordenadas normalizadas (0…1 sobre el ancho y alto del video)
  for(let i=0;i<n;i++)for(let c=0;c<4;c++){S[i*8+2*c]=(S[i*8+2*c]+.5)/w;S[i*8+2*c+1]=(S[i*8+2*c+1]+.5)/h;}
  let found=0;for(let i=0;i<n;i++)if(valid(i))found++;
  return{t:Float64Array.from(T),q:S,found:n?found/n:0};
}

/**
 * Analiza el video y devuelve dónde está la superficie en cada fotograma.
 * M: módulo de Mediabunny. keys: [{t, q:[{x,y}×4] normalizadas}]. onProgress(0…1). cancelled() → true para abortar.
 */
export async function analyze({M,file,keys,onProgress=()=>{},cancelled=()=>false}){
  const input=new M.Input({source:new M.BlobSource(file),formats:M.ALL_FORMATS});
  const track=await input.getPrimaryVideoTrack();
  if(!track||!(await track.canDecode()))throw new Error('nodec');
  const dw=track.displayWidth,dh=track.displayHeight,s=Math.min(1,640/Math.max(dw,dh));
  const w=Math.max(48,Math.round(dw*s)),h=Math.max(48,Math.round(dh*s));
  const sink=new M.CanvasSink(track,{width:w,height:h,fit:'fill',poolSize:2});
  const t0=await track.getFirstTimestamp(),t1=await track.computeDuration(),span=Math.max(1e-3,t1-t0);
  const stop=()=>{if(cancelled()){const e=new Error('cancel');e.name='Cancel';throw e;}};
  // Referencias: el fotograma exacto de cada momento marcado. Todas se expresan en la parametrización del
  // primero (keys[0]): así una corrección que solo arregla el seguimiento no deforma el resto del video,
  // y una que de verdad mueve las esquinas se reparte suavemente entre momentos marcados.
  const refs=[],keyMap=new Map(),keyT=new Set(),keyU=[],rand=rng(7);
  for(const k of keys){
    const wc=await sink.getCanvas(Math.max(t0,k.t));if(!wc)continue;
    const gray=grayOf(wc.canvas,w,h),quad=k.q.flatMap(p=>[p.x*w-.5,p.y*h-.5]),Hk=hFromQuad(quad);
    let ref=buildRef(gray,w,h,quad),C=null;
    if(refs.length){const H0k=relate(refs[0],ref,rand);if(H0k)C=hMul(hInv(Hk),H0k);}
    let U=[0,0,1,0,1,1,0,1];
    if(C){const Ci=hInv(C);U=[...hApply(Ci,0,0),...hApply(Ci,1,0),...hApply(Ci,1,1),...hApply(Ci,0,1)];
      if(U.every(Number.isFinite)&&Math.max(...U.map(Math.abs))<5)ref=buildRef(gray,w,h,hQuad(hMul(Hk,C)));else U=[0,0,1,0,1,1,0,1];}
    refs.push(ref);keyU.push({t:wc.timestamp,U});
    keyMap.set(tk(wc.timestamp),{H:ref.H,ref:refs.length-1});keyT.add(tk(wc.timestamp));
  }
  if(!refs.length)throw new Error('nokey');
  stop();
  // pasada hacia adelante (máx. ~30 fotogramas por segundo; el resto se interpola)
  const T=[],proc=new Set(),PT=[],fw=new Map(),bw=new Map();
  let tr=new Tracker(w,h,refs),last=-1e9;
  for await(const wc of sink.canvases()){
    stop();
    const ts=wc.timestamp,k=tk(ts);T.push(ts);
    if(!keyMap.has(k)&&ts-last<1/32-1e-4)continue;
    last=ts;proc.add(k);PT.push(ts);
    fw.set(k,tr.step(grayOf(wc.canvas,w,h),keyMap.get(k)));
    onProgress(.5*(ts-t0)/span);
  }
  // pasada hacia atrás, por tramos de ~2 s para no guardar todo el video en memoria
  tr=new Tracker(w,h,refs);
  let end=PT.length-1;
  while(end>=0){
    let start=end;while(start>0&&PT[end]-PT[start-1]<2)start--;
    const want=new Set(PT.slice(start,end+1).map(tk)),frames=[];
    for await(const wc of sink.canvases(PT[start],PT[end]+1e-4)){
      stop();const k=tk(wc.timestamp);
      if(want.has(k)&&!bw.has(k)&&!frames.some(f=>f.k===k))frames.push({k,g:grayOf(wc.canvas,w,h)});
    }
    for(let i=frames.length-1;i>=0;i--){const f=frames[i];stop();bw.set(f.k,tr.step(f.g,keyMap.get(f.k)));}
    onProgress(.5+.5*(t1-PT[start])/span);
    end=start-1;
  }
  T.sort((a,b)=>a-b);
  return merge(T,proc,fw,bw,keyT,w,h,keyU);
}
