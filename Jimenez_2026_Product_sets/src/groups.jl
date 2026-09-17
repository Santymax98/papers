# Hölder groups K(m,n,g,h), with τ^n=1, η^m=τ^g, ητη⁻¹=τ^h.
# Only integer arithmetic and Julia standard libraries are used.

"Prime factorization of a positive integer."
function factor_integer(n::Integer)
    n > 0 || throw(ArgumentError("n must be positive"))
    x = Int(n)
    out = Pair{Int,Int}[]
    p = 2
    while p <= x ÷ p
        if x % p == 0
            e = 0
            while x % p == 0
                x ÷= p
                e += 1
            end
            push!(out, p => e)
        end
        p = p == 2 ? 3 : p + 2
    end
    x > 1 && push!(out, x => 1)
    out
end

"Positive divisors in increasing order."
function divisors(n::Integer)
    n > 0 || throw(ArgumentError("n must be positive"))
    ds = Int[1]
    for (p, e) in factor_integer(n)
        previous = copy(ds)
        pe = 1
        for _ in 1:e
            pe *= p
            append!(ds, (d * pe for d in previous))
        end
    end
    sort!(ds)
end

function euler_phi(n::Integer)
    n > 0 || throw(ArgumentError("n must be positive"))
    value = Int(n)
    for (p, _) in factor_integer(n)
        value = value ÷ p * (p - 1)
    end
    value
end

function is_prime_integer(n::Integer)
    n >= 2 || return false
    n == 2 && return true
    iseven(n) && return false
    p = 3
    while p <= n ÷ p
        n % p == 0 && return false
        p += 2
    end
    true
end

"Multiplicative order of a modulo n."
function multiplicative_order(a::Integer, n::Integer)
    n > 1 || throw(ArgumentError("n must exceed 1"))
    gcd(a, n) == 1 || throw(ArgumentError("a must be a unit modulo n"))
    order = euler_phi(n)
    for (q, e) in factor_integer(order)
        for _ in 1:e
            candidate = order ÷ q
            powermod(mod(a, n), candidate, n) == 1 || break
            order = candidate
        end
    end
    order
end

struct MetacyclicGroup
    m::Int
    n::Int
    g::Int
    h::Int
    function MetacyclicGroup(m::Integer, n::Integer, g::Integer, h::Integer)
        m > 0 && n > 0 || throw(ArgumentError("m and n must be positive"))
        mm, nn = Int(m), Int(n)
        gg, hh = mod(Int(g), nn), mod(Int(h), nn)
        gcd(hh, nn) == 1 || throw(ArgumentError("h must be a unit modulo n"))
        mod(gg * (hh - 1), nn) == 0 ||
            throw(ArgumentError("the condition n | g(h-1) fails"))
        mod(powermod(hh, mm, nn) - 1, nn) == 0 ||
            throw(ArgumentError("the condition n | h^m-1 fails"))
        new(mm, nn, gg, hh)
    end
end

Base.show(io::IO, G::MetacyclicGroup) =
    print(io, "K($(G.m),$(G.n),$(G.g),$(G.h))")

Base.:(==)(G::MetacyclicGroup, H::MetacyclicGroup) =
    (G.m, G.n, G.g, G.h) == (H.m, H.n, H.g, H.h)
Base.hash(G::MetacyclicGroup, h::UInt) = hash((G.m, G.n, G.g, G.h), h)

"Element τ^j η^i in Yang normal form."
struct MetacyclicElement
    group::MetacyclicGroup
    i::Int
    j::Int
    function MetacyclicElement(G::MetacyclicGroup, i::Integer, j::Integer)
        new(G, mod(Int(i), G.m), mod(Int(j), G.n))
    end
end

Base.show(io::IO, x::MetacyclicElement) =
    x.i == 0 && x.j == 0 ? print(io, "e") :
    x.i == 0 ? print(io, "τ^", x.j) :
    x.j == 0 ? print(io, "η^", x.i) : print(io, "τ^", x.j, "η^", x.i)

Base.:(==)(x::MetacyclicElement, y::MetacyclicElement) =
    x.group == y.group && x.i == y.i && x.j == y.j
Base.hash(x::MetacyclicElement, h::UInt) = hash((x.group, x.i, x.j), h)

identity_element(G::MetacyclicGroup) = MetacyclicElement(G, 0, 0)

function Base.:*(x::MetacyclicElement, y::MetacyclicElement)
    x.group == y.group || throw(ArgumentError("elements belong to different groups"))
    G = x.group
    isum = x.i + y.i
    i = mod(isum, G.m)
    j = mod(x.j + powermod(G.h, x.i, G.n) * y.j + fld(isum, G.m) * G.g, G.n)
    MetacyclicElement(G, i, j)
end

function Base.inv(x::MetacyclicElement)
    G = x.group
    x.i == 0 && return MetacyclicElement(G, 0, -x.j)
    i = G.m - x.i
    j = -powermod(G.h, i, G.n) * (x.j + G.g)
    MetacyclicElement(G, i, j)
end

elements(G::MetacyclicGroup) =
    [MetacyclicElement(G, i, j) for i in 0:G.m-1 for j in 0:G.n-1]

element_index(x::MetacyclicElement) = x.i * x.group.n + x.j + 1

"Construct C_p ⋊ C_m with an action of exact order m."
function semidirect_group(p::Integer, m::Integer; h::Union{Nothing,Integer}=nothing)
    pp, mm = Int(p), Int(m)
    is_prime_integer(pp) || throw(ArgumentError("p must be prime"))
    mm > 0 && (pp - 1) % mm == 0 ||
        throw(ArgumentError("m must divide p-1"))

    action = if h === nothing
        found = findfirst(a -> multiplicative_order(a, pp) == mm, 1:pp-1)
        found === nothing && error("no action of order m was found modulo p")
        found
    else
        mod(Int(h), pp)
    end
    multiplicative_order(action, pp) == mm ||
        throw(ArgumentError("h does not have multiplicative order m modulo p"))
    MetacyclicGroup(mm, pp, 0, action)
end

struct YangTriple
    order::Int
    ell::Int
    beta::Int
    normal::Bool
end

function _sum_h_powers(G::MetacyclicGroup, k::Int, ell::Int)
    d = G.n ÷ ell
    d == 1 && return 0
    q = k ÷ ell
    step = G.m * ell ÷ k
    ratio = powermod(mod(G.h, d), step, d)
    power = mod(1, d)
    total = 0
    for _ in 1:q
        total = mod(total + power, d)
        power = mod(power * ratio, d)
    end
    total
end

function _is_normal_triple(G::MetacyclicGroup, k::Int, ell::Int, beta::Int)
    d = G.n ÷ ell
    first_condition = mod(beta * (G.h - 1), d) == 0
    second_condition = d == 1 || powermod(mod(G.h, d), G.m * ell ÷ k, d) == mod(1, d)
    first_condition && second_condition
end

"Yang triples parametrizing all subgroups of K(m,n,g,h)."
function yang_triples(G::MetacyclicGroup)
    triples = YangTriple[]
    for ell in divisors(G.n)
        d = G.n ÷ ell
        for k in divisors(G.m * ell)
            k >= ell && k % ell == 0 || continue
            S = _sum_h_powers(G, k, ell)
            for beta in 0:d-1
                mod(beta * S + G.g, d) == 0 || continue
                push!(triples, YangTriple(k, ell, beta,
                    _is_normal_triple(G, k, ell, beta)))
            end
        end
    end
    triples
end

function subgroup_elements(G::MetacyclicGroup, triple::YangTriple)
    k, ell, beta = triple.order, triple.ell, triple.beta
    q = k ÷ ell
    step = G.n ÷ ell
    u = MetacyclicElement(G, 0, beta) * MetacyclicElement(G, G.m * ell ÷ k, 0)
    H = Set{MetacyclicElement}()
    upower = identity_element(G)
    for a in 0:q-1
        a > 0 && (upower *= u)
        for s in 0:ell-1
            push!(H, upower * MetacyclicElement(G, 0, s * step))
        end
    end
    length(H) == k || error("invalid Yang subgroup construction")
    sort!(collect(H); by=element_index)
end

function subgroup_orders(G::MetacyclicGroup)
    triples = yang_triples(G)
    all_orders = sort!(unique([t.order for t in triples]))
    normal_orders = sort!(unique([t.order for t in triples if t.normal]))
    all_orders, normal_orders
end

key(G)=(G.m,G.n,G.g,G.h)

function subgroup_parameters(G,tr)
    ell=tr.ell; q=tr.order÷ell; D=G.n÷ell; t=G.m÷q
    S=sum(powermod(G.h,t*j,G.n) for j=0:q-1)
    residue=mod(G.g+tr.beta*S,G.n)
    @assert residue%D==0
    MetacyclicGroup(q,ell,residue÷D,powermod(G.h,t,ell))
end

function quotient_parameters(G,tr)
    @assert tr.normal
    q=tr.order÷tr.ell
    MetacyclicGroup(G.m÷q,G.n÷tr.ell,-tr.beta,G.h)
end

function cyclic_quotient(G,tr)
    tr.normal || return false
    D=G.n÷tr.ell; t=G.m÷(tr.order÷tr.ell)
    mod(G.h-1,D)==0 && gcd(gcd(D,t),tr.beta)==1
end


product_bound(r,s,d) = d * (cld(r,d)+cld(s,d)-1)
