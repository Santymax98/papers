# Cardinality implementation of the arithmetic trio region in the manuscript.
# U[r,s] = |G| - max{t : (r,s,t) in T_G}; zero coordinates are implicit.
# Noncritical vectors are already covered by the elementary construction.
include(joinpath(@__DIR__, "groups.jl"))

function dual_close!(U)
    N=size(U,1); changed=true; rounds=0
    while changed
        changed=false;rounds+=1
        for r=N:-1:1,s=N:-1:1
            old=U[r,s]; new=min(old,U[s,r])
            r<N && (new=min(new,U[r+1,s]))
            s<N && (new=min(new,U[r,s+1]))
            if new<old; U[r,s]=new;changed=true;end
        end
        for s=1:N,t=1:N
            r=N-U[s,t]
            if r>0 && N-t<U[r,s]
                U[r,s]=N-t;changed=true
            end
        end
    end
    rounds
end

const devos_cache = Dict{NTuple{4,Int},Matrix{Int}}()
const boundary_cache = Dict{NTuple{4,Int},Vector{NTuple{3,Int}}}()
const dihedral_cache = Dict{NTuple{4,Int},Int}()
const PERMS3 = ((1,2,3),(1,3,2),(2,1,3),(2,3,1),(3,1,2),(3,2,1))
const OCT_FACES = [(2*u+v+1, 5+2*v+w, 9+2*w+u)
                   for u=0:1 for v=0:1 for w=0:1]

function mask_counts(kind)
    grey = kind==0 ? Set{Int}() : kind==1 ? Set(OCT_FACES[1]) :
           union(Set(OCT_FACES[1]),Set(OCT_FACES[2]))
    allowed = kind==0 ? Int[] : kind==1 ? [1] : [1,2]
    remaining = [i for i=1:12 if !(i in grey)]
    out=Set{NTuple{3,Int}}()
    for mask=0:(1<<length(remaining))-1
        full=Set(remaining[j] for j in eachindex(remaining) if (mask>>(j-1))&1==1)
        nz=union(grey,full)
        any(!(i in allowed) && all(e in nz for e in face)
            for (i,face) in enumerate(OCT_FACES)) && continue
        push!(out,ntuple(j->count(e->e in full,4*j-3:4*j),3))
    end
    # Coordinatewise maximal counts suffice and greatly reduce arithmetic work.
    [x for x in out if !any(x!=y && all(x[j]<=y[j] for j=1:3) for y in out)]
end
const OCT_COUNTS = [mask_counts(k) for k=0:2]

function seed_trio!(U,r,s,t)
    N=size(U,1)
    (1<=r<=N && 1<=s<=N && 1<=t<=N) || return
    v=(r,s,t)
    for p in PERMS3
        a,b,c=v[p[1]],v[p[2]],v[p[3]]
        U[a,b]=min(U[a,b],N-c)
    end
end

function critical_rows(U)
    D=size(U,1)
    [(a,b,D-U[a,b]) for a=1:D for b=1:D
        if U[a,b]<D && a+b>U[a,b]]
end

function dihedral_degree(G)
    get!(dihedral_cache,key(G)) do
        N=G.m*G.n
        (iseven(N) && N>=6) || return 0
        n=N÷2; E=elements(G); e=identity_element(G)
        for r in E
            R=Set([e]); z=r
            while !(z in R)
                push!(R,z); z=z*r
            end
            length(R)==n || continue
            for f in E
                if !(f in R) && f*f==e && f*r*f==inv(r)
                    return n
                end
            end
        end
        0
    end
end

function boundary_totals(H)
    haskey(boundary_cache,key(H)) && return boundary_cache[key(H)]
    D=H.m*H.n; out=Set{NTuple{3,Int}}()
    function add(f,proper)
        v=ntuple(j->D*f[j]+proper[j],3)
        sum(v)>6D || return # only critical boundaries are needed
        for p in PERMS3
            push!(out,(v[p[1]],v[p[2]],v[p[3]]))
        end
    end
    for f in OCT_COUNTS[1];add(f,(0,0,0));end
    for (a,b,c) in critical_rows(devos_table(H)), f in OCT_COUNTS[2]
        add(f,(a,b,c))
    end
    trs=[tr for tr in yang_triples(H) if tr.order<D]
    sets=[Set(subgroup_elements(H,tr)) for tr in trs]
    for (i,trK) in enumerate(trs), (j,trL) in enumerate(trs)
        k=trK.order; ell=trL.order
        # Type 2B: two pure beats with a common intersection component.
        a=length(intersect(sets[i],sets[j]))
        for u=1:D÷k-1, v=1:D÷ell-1, f in OCT_COUNTS[3]
            add(f,(a,u*k+v*ell,2D-u*k-v*ell))
        end
        # Type 2A: an impure beat and a containing pure beat.
        issubset(sets[i],sets[j]) || continue
        for (a,b,c) in critical_rows(devos_table(subgroup_parameters(H,trK)))
            for u=0:D÷k-1, v=1:D÷ell-1, f in OCT_COUNTS[3]
                add(f,(a,u*k+b+v*ell,2D-(u+1)*k+c-v*ell))
            end
        end
    end
    # Downward domination is safe for subsequent use of just the three totals.
    vals=collect(out)
    ans=[x for x in vals if !any(x!=y && all(x[j]<=y[j] for j=1:3) for y in vals)]
    boundary_cache[key(H)]=ans
end

function coset_seeds!(U,G,tr)
    N=size(U,1);D=tr.order
    pairs=[(a,b) for (a,b) in ((2,2),(2,3)) if b*D<=N &&
        U[a*D,b*D]>minimum(product_bound(a*D,b*D,d) for d in divisors(N))]
    isempty(pairs) && return
    E=elements(G); HS=subgroup_elements(G,tr)
    # Only representatives of at most three cosets are chosen, never arbitrary subsets.
    L=Vector{typeof(E[1])}(); R=Vector{typeof(E[1])}()
    seenL=Set{typeof(E[1])}(); seenR=Set{typeof(E[1])}()
    for x in E
        if !(x in seenL);push!(L,x);union!(seenL,(x*z for z in HS));end
        if !(x in seenR);push!(R,x);union!(seenR,(z*x for z in HS));end
    end
    blocks=[BitSet(element_index(x*z*y) for z in HS) for x in L,y in R]
    q=length(L)
    for (a,b) in pairs
        for ix=2:q, jy=2:q
            base=union(blocks[1,1],blocks[ix,1],blocks[1,jy],blocks[ix,jy])
            if b==2
                seed_trio!(U,2D,2D,N-length(base))
            else
                for jz=jy+1:q
                    P=union(base,blocks[1,jz],blocks[ix,jz])
                    seed_trio!(U,2D,3D,N-length(P))
                end
            end
        end
    end
end

function devos_table(G)
    haskey(devos_cache,key(G)) && return devos_cache[key(G)]
    N=G.m*G.n
    U=[min(N,r+s-1) for r=1:N,s=1:N]
    for tr in yang_triples(G)
        D=tr.order;D==N && continue
        H=subgroup_parameters(G,tr); V=devos_table(H)
        for r=1:D,s=1:N
            b,rem=divrem(s-1,D)
            U[r,s]=min(U[r,s],b*D+V[r,rem+1])
            U[s,r]=min(U[s,r],U[r,s])
        end
        if tr.normal
            Q=quotient_parameters(G,tr)
            if D>1
                W=devos_table(Q)
                for r=1:N,s=1:N
                    U[r,s]=min(U[r,s],D*W[cld(r,D),cld(s,D)])
                end
            end
            if cyclic_quotient(G,tr)
                for r=1:N,s=1:N
                    a,ra=divrem(r-1,D);b,rb=divrem(s-1,D)
                    U[r,s]=min(U[r,s],min(N,(a+b)*D+V[ra+1,rb+1]))
                end
            end
            n=dihedral_degree(Q)
            if n>0
                for k=0:n-2,l=0:n-k-2
                    seed_trio!(U,2*(k+1)*D,2*(l+1)*D,2*(n-k-l-1)*D)
                end
            end
            if n>=4
                for (x,y,z) in boundary_totals(H), k=1:n-3,l=1:n-k-2
                    seed_trio!(U,2*(k-1)*D+x,2*(l-1)*D+y,2*(n-k-l-1)*D+z)
                end
            end
        end
        coset_seeds!(U,G,tr)
    end
    dual_close!(U)
    ds=divisors(N)
    @assert all(U[r,s]>=minimum(product_bound(r,s,d) for d in ds) for r=1:N,s=1:N)
    devos_cache[key(G)]=U
end
