import numpy as np, itertools
q=11
# quadratic character
sq={ (x*x)%q for x in range(1,q) }
chi=np.array([0]+[1 if x in sq else -1 for x in range(1,q)])
xs=np.arange(q)
# P monic of degree 5 with P(0)=0 : coeffs a1..a4
sigs=set()
pw={d: (xs[:,None]**d % q) for d in range(7)}
for a in itertools.product(range(q),repeat=4):
    P=(pw[5][:,0]+a[3]*pw[4][:,0]+a[2]*pw[3][:,0]+a[1]*pw[2][:,0]+a[0]*xs)%q
    # N_c = q + sum_x chi(P(x)+c)
    N=np.array([q+int(chi[(P+c)%q].sum()) for c in range(q)])
    sigs.add(tuple(sorted(N.tolist(),reverse=True)))
# degree 6 monic, P(0)=0: a1..a5
for a in itertools.product(range(q),repeat=5):
    P=(pw[6][:,0]+a[4]*pw[5][:,0]+a[3]*pw[4][:,0]+a[2]*pw[3][:,0]+a[1]*pw[2][:,0]+a[0]*xs)%q
    N=np.array([q+int(chi[(P+c)%q].sum()) for c in range(q)])
    sigs.add(tuple(sorted(N.tolist(),reverse=True)))
sigs=sorted(sigs)
print('distinct signatures:',len(sigs))
np.save('f11_sigs.npy', np.array(sigs))
