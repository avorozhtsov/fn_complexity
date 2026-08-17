import numpy as np, itertools, math
S=np.load("f11_sigs.npy"); q=11
S=S[(S>0).all(axis=1)]   # signatures must be positive integers; drop empty fibers
n=len(S); L=np.log(S.astype(float))                    # n x q
# beta grid: brief B says horizon must reach ~360q ; include 0 and infinity separately
b=np.concatenate([np.linspace(0,2,600)[1:], np.geomspace(2,360*q,3400)])
# logZ(beta) = logsumexp(beta * log N_c)
Zl=np.empty((n,len(b)))
for i in range(n):
    z=np.outer(b,L[i])                                  # len(b) x q
    m=z.max(axis=1); Zl[i]=m+np.log(np.exp(z-m[:,None]).sum(axis=1))
U=np.log(Zl)                                            # u = log log Z
u0=np.log(np.log(float(q)))                             # beta=0: Z=q fibers
uinf=np.log(L.max(axis=1))                              # beta=inf: u -> log(beta*log max)+..., difference -> log(logmaxa/logmaxb)
# A(a,b) = mid over beta of (u_a - u_b); endpoints: beta=0 gives 0, beta=inf gives uinf_a-uinf_b
A=np.zeros((n,n))
for i in range(n):
    D=U[i][None,:]-U                                    # n x len(b)
    lo=np.minimum(D.min(axis=1), np.minimum(0.0, uinf[i]-uinf))
    hi=np.maximum(D.max(axis=1), np.maximum(0.0, uinf[i]-uinf))
    A[i]=0.5*(lo+hi)
A=(A-A.T)/2
psi=-A.mean(axis=1); G=psi[None,:]-psi[:,None]
nA=np.linalg.norm(A)
print(f'n={n}  ||A||={nA:.4f}  ||grad||/||A||={np.linalg.norm(G)/nA:.4f}  ||curl||/||A||={np.linalg.norm(A-G)/nA:.4f}')
tol=1e-9
P=A>tol
cyc=sum(1 for i,j,k in itertools.combinations(range(n),3)
        if (P[i,j] and P[j,k] and P[k,i]) or (P[j,i] and P[k,j] and P[i,k]))
print('strict 3-cycles:',cyc,'of',math.comb(n,3))
np.save('f11_A.npy',A); np.save('f11_psi.npy',psi)
# what is psi?  regress against M = max_c(-a_c) = maxN - q, m2, mult of max
maxN=S.max(axis=1).astype(float); M=maxN-q
a_c=q-S.astype(float); m2=(a_c**2).sum(axis=1)/q**2
m3=(a_c**3).sum(axis=1)/q**3; m4=(a_c**4).sum(axis=1)/q**4
mult=(S==S.max(axis=1)[:,None]).sum(axis=1).astype(float)
phi=np.log(float(q))*np.log(maxN)
def reg(cols,names):
    X=np.column_stack([np.ones(n)]+cols); beta,*_=np.linalg.lstsq(X,psi,rcond=None)
    r=psi-X@beta; R2=1-r.var()/psi.var()
    print(f'  R2={R2:.6f}  {names}')
print('psi_opt regressions:')
reg([M],'M')
reg([np.log(maxN)],'log maxN')
reg([0.5*np.log(phi)],'(1/2)log phi')
reg([M,m2],'M,m2')
reg([M,m2,np.log(mult)],'M,m2,log mult')
reg([M,m2,m3,m4,np.log(mult)],'M,m2,m3,m4,log mult')
