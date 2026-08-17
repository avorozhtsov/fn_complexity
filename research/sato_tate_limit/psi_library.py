import numpy as np, itertools
NB=4000
def hist(alpha,w,lo,hi):
    h,e=np.histogram(alpha,bins=NB,range=(lo,hi),weights=w)
    ctr=0.5*(e[1:]+e[:-1]); m=h>0
    return ctr[m], h[m]/h.sum()
class M:
    def __init__(s,a,w,amax): s.a,s.w,s.amax=a,w/w.sum(),amax
    def K(s,t):
        z=t*s.a; m=z.max(); return np.log(np.sum(s.w*np.exp(z-m)))+m
th=np.linspace(0,np.pi,20001); c=2*np.cos(th)
sc=(2/np.pi)*np.sin(th)**2; un=np.ones_like(th)/len(th)
def conv(a1,w1,a2,w2,amax):
    A,B=np.meshgrid(a1,a2,indexing='ij')
    a,w=hist((A+B).ravel(),np.outer(w1,w2).ravel(),-amax,amax)
    return M(a,w,amax)
# base measures binned first
su2=hist(c,sc,-2,2); u1=hist(c,un,-2,2)
cm_a=np.concatenate([np.zeros(10001),c]); cm_w=np.concatenate([np.full(10001,.5/10001),np.full(len(th),.5/len(th))])
cm=hist(cm_a,cm_w,-2,2)
lib={}
g=np.linspace(0,np.pi,600);T1,T2=np.meshgrid(g,g,indexing='ij')
a,w=hist((2*np.cos(T1)+2*np.cos(T2)).ravel(),((np.cos(T1)-np.cos(T2))**2*np.sin(T1)**2*np.sin(T2)**2).ravel(),-4,4)
lib['USp4']=M(a,w,4.0)
lib['SU2xSU2']=conv(*su2,*su2,4.0)
lib['U1xU1']=conv(*u1,*u1,4.0)
lib['SU2xU1']=conv(*su2,*u1,4.0)
lib['CMxCM']=conv(*cm,*cm,4.0)
lib['SU2xCM']=conv(*su2,*cm,4.0)
lib['U1xCM']=conv(*u1,*cm,4.0)
lib['SU2_mult2']=M(2*su2[0],su2[1],4.0)
lib['U1_mult2']=M(2*u1[0],u1[1],4.0)
lib['CM_mult2']=M(2*cm[0],cm[1],4.0)
names=list(lib)
tau=np.concatenate([np.linspace(1e-3,2,300),np.linspace(2.01,20,500),np.linspace(20.1,400,500)])
P={n:(np.array([lib[n].K(t)/t for t in tau]),lib[n].amax) for n in names}
def mid(i,j):
    a,ai=P[i]; b,bj=P[j]
    D=np.concatenate([[0.0],a-b,[ai-bj]]); return 0.5*(D.max()+D.min())
print('             '+''.join(f'{n:>11s}' for n in names))
for i in names: print(f'{i:12s}'+''.join(f'{mid(i,j):11.6f}' for j in names))
cyc=[]
for t in itertools.permutations(names,3):
    x,y,z=t
    if x==min(t) and mid(x,y)<0 and mid(y,z)<0 and mid(z,x)<0: cyc.append((t,mid(x,y),mid(y,z),mid(z,x)))
print('CYCLES:',len(cyc))
for t,m1,m2,m3 in cyc: print('   ',t,f'{m1:.4e} {m2:.4e} {m3:.4e}')

print()
print('=== crossing test: does Psi_mu - Psi_nu change sign on (0,inf)? ===')
ncross=0
for i,j in itertools.combinations(names,2):
    a,ai=P[i]; b,bj=P[j]
    D=np.concatenate([[0.0],a-b,[ai-bj]])
    Dn=D[1:-1]
    if Dn.max()>1e-9 and Dn.min()<-1e-9:
        ncross+=1
        print(f'  CROSS {i:11s} vs {j:11s} range=[{Dn.min():+.5f},{Dn.max():+.5f}] mid={0.5*(D.max()+D.min()):+.5f}')
print('pairs with crossing Psi:',ncross,'of',len(names)*(len(names)-1)//2)
print()
print('=== is mid a potential difference?  fit psi and report residual ===')
import numpy as _np
n=len(names); Kmat=_np.zeros((n,n))
for x in range(n):
    for y in range(n): Kmat[x,y]=mid(names[x],names[y])
psi=-Kmat.mean(axis=1)
G=psi[None,:]-psi[:,None]
print('  ||K||=%.6f  ||K-grad||/||K|| = %.3e' % (_np.linalg.norm(Kmat), _np.linalg.norm(Kmat-G)/_np.linalg.norm(Kmat)))
print('  psi values:', {names[k]: round(float(psi[k]),6) for k in range(n)})
