include(joinpath(@__DIR__, "..", "src", "recursion.jl"))
using Test

power(x,k)=foldl(*,fill(x,k);init=identity_element(x.group))
productset(A,B)=Set(a*b for a in A for b in B)

@testset "Section 7 examples" begin
    println("p  m  r  D_kappa  mu  N_kappa  witness")
    for (p,m,r,expected) in ((17,4,8,(12,14,15)),(31,6,12,(18,22,23)),
                             (37,6,18,(30,34,35)))
        G=semidirect_group(p,m); U=devos_table(G)
        ds,nds=subgroup_orders(G)
        lower=minimum(product_bound(r,r,d) for d in divisors(p*m))
        upper=minimum(product_bound(r,r,d) for d in nds)
        @test ds==divisors(p*m)
        @test (lower,U[r,r],upper)==expected
        I=[MetacyclicElement(G,0,j) for j=0:r÷2-1]
        F=[identity_element(G),MetacyclicElement(G,m÷2,0)]
        A=productset(I,F); B=productset(F,I); AB=productset(A,B)
        @test length(A)==length(B)==r
        @test length(AB)==U[r,r]
        println("$p  $m  $r  $lower  $(U[r,r])  $upper  $(length(AB))")
    end

    G=MetacyclicGroup(9,57,19,4); U=devos_table(G)
    eta=MetacyclicElement(G,1,0); tau=MetacyclicElement(G,0,1)
    N=Set(x for x in elements(G) if x.i%3==0)
    H=Set(power(eta,3*j) for j=0:8)
    @test length(Set(power(eta,j) for j=0:26))==27 && power(eta,27)==identity_element(G)
    @test length(H)==9 && length(N)==171
    @test all(power(power(tau,b)*eta,9)==power(tau,19)!=identity_element(G) for b=0:56)
    A=union(N,Set(eta*x for x in H)); AB=productset(A,A)
    @test length(A)==180
    @test AB==union(N,Set(eta*x for x in N),Set(eta*eta*x for x in H))
    lower=minimum(product_bound(180,180,d) for d in divisors(513))
    @test lower==U[180,180]==length(AB)==351
    println("K(9,57,19,4): r=s=180, D_kappa=$lower, mu=$(U[180,180]), witness=$(length(AB))")
end
