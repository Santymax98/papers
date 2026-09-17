"""Independent integer-coordinate tests of the dihedral boundary lift.
No Julia code, profile bound, or claimed exact formula is used.
"""
from itertools import product
import json

FACES=[(2*u+v,4+2*v+w,8+2*w+u) for u,v,w in product(range(2),repeat=3)]

def masks(kind):
    proper=set() if kind==0 else set(FACES[0]) if kind==1 else set(FACES[0]+FACES[1])
    allowed=set(range(kind+0)) if kind<2 else {0,1}
    rem=[e for e in range(12) if e not in proper]
    for bits in product(range(2),repeat=len(rem)):
        full={e for e,b in zip(rem,bits) if b}
        if any(i not in allowed and set(face)<=proper|full for i,face in enumerate(FACES)):
            continue
        yield full

def audit(D,kind,split,q=4,k=1,l=1):
    n=q*D;g=0 if split else n//2;h=n-1
    def mul(a,b):
        i,x=a;j,y=b
        return ((i+j)%2,(pow(h,j,n)*x+y+g*((i+j)//2))%n)
    def inv(a):
        i,x=a
        return (0,(-x)%n) if i==0 else (1,(-h*x-g)%n)
    E={(i,x) for i in range(2) for x in range(n)};one=(0,0)
    H={(0,q*z) for z in range(D)}
    # Calculate the quotient prechord independently of the Julia recursion.
    def qm(a,b):
        i,x=a;j,y=b
        return ((i+j)%2,((-1)**j*x+y)%q)
    def qi(a):
        i,x=a
        return (i,(-(-1)**i*x)%q)
    QE={(i,x) for i in range(2) for x in range(q)}
    r=(0,1);f=(1,0);qone=(0,0)
    Ap={qm(x,y) for x in [(0,i) for i in range(k+1)] for y in (qone,f)}
    Bp={qm(x,y) for x in (qone,f) for y in [(0,j) for j in range(l+1)]}
    XP={qm(x,y) for x in Ap for y in Bp}
    endsC={qone,qm(f,(0,-k)),(0,(-k-l)%q),qm(f,(0,l))}
    Cp=(QE-{qi(x) for x in XP})|endsC
    if kind==0: labels={}
    elif kind==1: labels=dict(zip(FACES[0],({one},{one},H-{one})))
    elif kind==2:
        K={(0,n//2*z) for z in range(2)}
        L={(0,n//4*z) for z in range(4)}
        P={one};Q={one};R=H-{one};S=L;T=H-L
        labels=dict(zip((0,4,8,5,10),(P,Q,R,S,T)))
    else:
        K={(0,n//2*z) for z in range(2)}
        L={(0,n//3*z) for z in range(3)}
        P=K&L;Q=K;R=H-K;S=L;T=H-L
        labels=dict(zip((0,4,8,5,10),(P,Q,R,S,T)))
    qtype=min(kind,2);count=0
    for lift_variant in (0,1):
        # Change three representatives by kernel elements when variant=1.
        X=[one,mul(mul(f,(0,(-k)%n)),(0,q*lift_variant))]
        Y=[one,mul(f,(0,2*q*lift_variant))]
        Z=[one,mul(mul(f,(0,l)),(0,3*q*lift_variant))]
        vertices=[(X[u],Y[v]) for u,v in product(range(2),repeat=2)]
        vertices += [(Y[v],Z[w]) for v,w in product(range(2),repeat=2)]
        vertices += [(Z[w],X[u]) for w,u in product(range(2),repeat=2)]
        qe=[(mul(inv(x),y)[0],mul(inv(x),y)[1]%q) for x,y in vertices]
        assert len(set(qe[:4]))==len(set(qe[4:8]))==len(set(qe[8:]))==4
        hulls=[Ap,Bp,Cp]
        interior=[{e for e in E if (e[0],e[1]%q) in hulls[j]-set(qe[4*j:4*j+4])}
                  for j in range(3)]
        for full in masks(qtype):
            labs=[labels.get(e,H if e in full else set()) for e in range(12)]
            sets=[set(s) for s in interior]
            for e,(x,y) in enumerate(vertices):
                sets[e//4].update(mul(mul(inv(x),z),y) for z in labs[e])
            A,B,C=sets
            AB={mul(a,b) for a in A for b in B}
            assert AB.isdisjoint({inv(c) for c in C})
            totals=[sum(len(labs[e]) for e in range(4*j,4*j+4)) for j in range(3)]
            assert [len(s) for s in sets]==[2*(k-1)*D+totals[0],2*(l-1)*D+totals[1],2*(q-k-l-1)*D+totals[2]]
            count+=1
    return dict(kernel_order=D,quotient_order=2*q,k=k,l=l,boundary_type=('0','1','2A','2B')[kind],split=split,
                verified_lifts=count)

# A second extension has a NONABELIAN normal subgroup H ~= D_13 and quotient D_5.
# G=K(4,65,0,44). Its quotient reflections lift to elements of order four.
def audit_nonabelian_kernel():
    m,n=4,65
    h=pow(44,-1,n) # η^i τ^x coordinates: inverse of the Yang action 44.
    def mul(a,b):
        i,x=a;j,y=b
        return ((i+j)%m,(pow(h,j,n)*x+y)%n)
    def inv(a):
        i,x=a;j=(-i)%m
        return (j,(-pow(h,j,n)*x)%n)
    E={(i,x) for i in range(m) for x in range(n)};one=(0,0)
    H={(i,5*z) for i in (0,2) for z in range(13)}
    assert any(mul(x,y)!=mul(y,x) for x in H for y in H)
    K={one,(2,0)};r=(0,1);f=(1,0)
    L={mul(mul(r,z),inv(r)) for z in K}
    assert K!=L and K&L=={one}
    def qm(a,b):
        i,x=a;j,y=b
        return ((i+j)%2,((-1)**j*x+y)%5)
    def qi(a):
        i,x=a;return (i,(-(-1)**i*x)%5)
    def pi(a):return (a[0]%2,a[1]%5)
    qone=(0,0);qr=(0,1);qf=(1,0)
    Ap={qm(x,y) for x in (qone,qr) for y in (qone,qf)}
    Bp={qm(x,y) for x in (qone,qf) for y in (qone,qr)}
    XP={qm(x,y) for x in Ap for y in Bp}
    Cp=({(i,x) for i in (0,1) for x in range(5)}-{qi(x) for x in XP})|{qone,qm(qf,qi(qr)),(0,3),qm(qf,qr)}
    X=[one,mul(f,inv(r))];Y=[one,f];Z=[one,mul(f,r)]
    vertices=[(X[u],Y[v]) for u,v in product(range(2),repeat=2)]
    vertices += [(Y[v],Z[w]) for v,w in product(range(2),repeat=2)]
    vertices += [(Z[w],X[u]) for w,u in product(range(2),repeat=2)]
    qe=[pi(mul(inv(x),y)) for x,y in vertices]
    interior=[{e for e in E if pi(e) in hull-set(qe[4*j:4*j+4])}
              for j,hull in enumerate((Ap,Bp,Cp))]
    total=0
    for kind in (1,2,3):
        if kind==1:labels=dict(zip(FACES[0],({one},{one},H-{one})))
        elif kind==2:labels=dict(zip((0,4,8,5,10),({one},{one},H-{one},K,H-K)))
        else:labels=dict(zip((0,4,8,5,10),(K&L,K,H-K,L,H-L)))
        for full in masks(min(kind,2)):
            labs=[labels.get(e,H if e in full else set()) for e in range(12)]
            sets=[set(s) for s in interior]
            for e,(x,y) in enumerate(vertices):
                sets[e//4].update(mul(mul(inv(x),z),y) for z in labs[e])
            A,B,C=sets;AB={mul(a,b) for a in A for b in B}
            assert AB.isdisjoint({inv(c) for c in C})
            totals=[sum(len(labs[e]) for e in range(4*j,4*j+4)) for j in range(3)]
            assert [len(s) for s in sets]==[totals[0],totals[1],104+totals[2]]
            total+=1
    return total

if __name__=='__main__':
    results=[audit(D,kind,split) for D,kind in ((4,0),(4,1),(8,2),(6,3))
             for split in (True,False)]
    results += [audit(4,kind,False,q,k,l) for q,k,l in ((5,1,2),(5,2,1),(6,2,2)) for kind in (0,1)]
    print(json.dumps(results,indent=2))
    print('total_lifts',sum(x['verified_lifts'] for x in results))
    extra=audit_nonabelian_kernel()
    print('nonabelian_kernel_lifts',extra)
    print('grand_total_lifts',sum(x['verified_lifts'] for x in results)+extra)
