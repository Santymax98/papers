include(joinpath(@__DIR__, "..", "src", "recursion.jl"))
using Test

# Independent exhaustive pair enumeration: normalize 1 into both input sets.
# Left translation of A and right translation of B preserve |AB|.
function brute_table(G)
    E=elements(G);N=length(E);N<=12 || error("exhaustion restricted to order 12")
    mult=[element_index(x*y) for x in E,y in E]
    out=fill(N,N,N);full=1<<N
    prod=zeros(UInt64,full);left=zeros(UInt64,N)
    for B=1:2:full-1
        sb=count_ones(B)
        for x=1:N
            z=UInt64(0)
            for y=1:N
                (B>>(y-1))&1==1 && (z |= UInt64(1)<<(mult[x,y]-1))
            end
            left[x]=z
        end
        prod[1]=0
        for A=1:full-1
            bit=trailing_zeros(A)
            prod[A+1]=prod[(A&(A-1))+1]|left[bit+1]
            if isodd(A)
                sa=count_ones(A);out[sa,sb]=min(out[sa,sb],count_ones(prod[A+1]))
            end
        end
    end
    out
end

@testset "all valid Holder presentations through order 12" begin
    presentations=0;distinct_tables=Dict{Any,Matrix{Int}}();cells=0
    for N=1:12,m in divisors(N)
        n=N÷m
        for h=0:n-1,g=0:n-1
            gcd(h,n)==1 && mod(powermod(h,m,n)-1,n)==0 && mod(g*(h-1),n)==0 || continue
            G=MetacyclicGroup(m,n,g,h);E=elements(G)
            sig=Tuple(element_index(x*y) for x in E for y in E)
            one=identity_element(G)
            @test all(x*inv(x)==one && inv(x)*x==one for x in E)
            @test all((x*y)*z==x*(y*z) for x in E,y in E,z in E)
            tau=MetacyclicElement(G,0,1); eta=MetacyclicElement(G,1,0)
            power(x,k)=foldl(*,fill(x,k);init=one)
            if m>1
                @test eta*tau*inv(eta)==power(tau,h)
                @test power(eta,m)==power(tau,g)
            end
            for tr in yang_triples(G)
                H=Set(subgroup_elements(G,tr))
                @test length(H)==tr.order && all(x*y in H for x in H,y in H)
                @test tr.normal == all(Set(x*y*inv(x) for y in H)==H for x in E)
            end
            exact=get!(distinct_tables,sig) do;brute_table(G);end
            U=devos_table(G)
            @test U==exact
            presentations+=1;cells+=N*N
        end
    end
    println((small_presentations=presentations,distinct_multiplication_tables=length(distinct_tables),ordered_cells=cells))
end

@testset "full Berchenko-Kogan tables" begin
    for (p,q) in ((7,3),(11,5),(13,3))
        G=semidirect_group(p,q);U=devos_table(G); ds,nds=subgroup_orders(G)
        for r=1:p*q,s=r:p*q
            lower=minimum(product_bound(r,s,d) for d in ds)
            upper=minimum(product_bound(r,s,d) for d in nds)
            expected=r>q && s>q && cld(r,q)+cld(s,q)<p ? upper : lower
            @test U[r,s]==expected
        end
    end
end

@testset "faithful complement formula in its proved nonsaturated domain" begin
    for (p,m) in ((13,4),(19,6),(17,8),(19,9),(31,10),(13,12))
        G=semidirect_group(p,m);U=devos_table(G);ds=divisors(p*m)
        for r=1:p,s=r:p-r+1
            expected=min(r,s)<=m ? minimum(product_bound(r,s,d) for d in ds) :
                r+s-(iseven(m)&&iseven(r)&&iseven(s) ? 2 : 1)
            @test U[r,s]==expected
        end
        println((faithful=key(G),agrees=true));flush(stdout)
    end
end

@testset "complete split and nonsplit dihedral-extension tables" begin
    for G in (MetacyclicGroup(2,12,0,7),MetacyclicGroup(2,16,0,15),
              MetacyclicGroup(2,16,8,15),MetacyclicGroup(2,9,0,8))
        U=devos_table(G);N=G.m*G.n;ds=divisors(N)
        @test U==[minimum(product_bound(r,s,d) for d in ds) for r=1:N,s=1:N]
        println((extension=key(G),kappa=true));flush(stdout)
    end
end

