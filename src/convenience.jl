function benchpar(;β = 0.98,χ = 0.001,α=1.0,γ=2.0,κ=0.048,d=0.2,
        lc=0.0,lb=0.0,ls=0.0,
        nx = 41,ns=2,nh=2,ne=1,np=3,nν=3,nϵ=1,
            ms = 0.0,mb=0.0,welfare = true,
            contract_type = "convex",)

    ## Initialize housing and ownership grids
    if ns == 1
        sgrd = [0.0] # Only renting
    elseif ns == 2
        sgrd = [0.0, 1.0] # Rent or own
    else
        sgrd = [0.0; 1.0; collect(range(0.5,length=ns-2,step=0.5/(ns-2)))]
        # Note that the grid is NOT sorted. 0 (first) = rent, 1 (second) = own, rest is partial ownershi
    end

    if nh > 1
        hgrd = range(start=1.0,stop=2.0,length=nh)
    else
        hgrd = [1.0]
    end
    shgrd = fill(true,(ns,nh))

    ## Initializing education grids
    _fa,_hdr = readdlm(joinpath(@__DIR__,"../data/modelinputs/income.csv"),',';header=true)
    if ne ==1
        σνs = sqrt.([0.012])
        egrd = ["All"]

        fa = _fa[:,4]
        @assert _hdr[:,4] == egrd

        CDFe_ini = range(1.0,stop=1.0,length=ne)
    else
        if ne == 3
            egrd = ["<HS", "HS", "Coll."]
            σνs = sqrt.([0.013, 0.011, 0.011])

            fa = _fa[:,1:3]
            @assert _hdr[1:3] == egrd

            CDFe_ini = [0.22, 0.22+0.53, 1.0] # Table I, column 1, FGG (2017 JF)
        end
    end

    ϵgrd = fill(1.0,(nϵ,76))
    πϵ = fill(1.0/nϵ,(nϵ,76))

    if np == 1
        pgrd = [150.0]
        πp = fill(1.0,1,1)
    else
        pgrd = [0.7,1.0,1.45]*150
        πp = [0.95 0.05 0.00; 0.01 0.98 0.01; 0.00 0.12 0.88]
    end
    ## Initializing initial x distribution
    nxini = 21
    dist = LogNormal(log(1.0),2.5) # Try also Normal(), Pareto()
    xgrd_ini,_,CDFx_ini = equiprob_discr(dist,nxini)

    CDFν_ini = fill(0.0,nν,ne)
    πν =  fill(0.0,100,nν,nν,ne)


    if welfare == true
        rentsupport = 10.0 # On top of rent
        sellersupport = 10.0 # On top of rent after selling
    else
        rentsupport = -Inf64 # On top of rent
        sellersupport = -Inf64 # On top of rent after mincost_selling_to_rent
    end

    par = (γ=γ,
        β = β,
        χ = χ,
        η = 0.25,
        κ = κ,
        d = d,
        b2y = 5.0,
        δ = 0.025,
        rf = 0.022,
        rm = 0.01,
        ρ = 0.95,
        α = α,

        ms = ms,
        mb = mb,

        ls = ls, #Fixed cost of selling shared ownership (s∈(0,1))
        lc = lc, #Fixed cost of changing ownershare (s∈(0,1)!=s' & h==h')
        lb = lb, #Fixed cost of buying shared ownership (s'∈(0,1))

        Ts = 25,
        Tr = 67,
        Te = 100,

        xgrd = range(start=-8.0,stop=1000.0,length=nx), # Wealth
        hgrd = hgrd, # Contains all house sizes
        sgrd = sgrd, # 0=rent, 1=own.
        shgrd = shgrd, # Grid for housing×ownership sets (row is tenure (s), column is house quality (h))
        egrd = egrd, # Discrete education grid
        νgrd = fill(1.0,nν), # Persistent Income Shock (productivity)
        ϵgrd = ϵgrd,
        pgrd = pgrd,
        agegrd = range(25,stop=100,step=1),
        xgrd_ini = xgrd_ini,

        πν = πν,
        πp = πp,
        πϵ = πϵ,
        CDFν_ini = CDFν_ini,
        CDFx_ini = CDFx_ini,
        CDFe_ini = CDFe_ini,

        fa = fa,

        nϵ = nϵ,

        ## Simulation
        N = 6000,
        seed = 221956, # Serial Number

        ## Various parameters
        minrentalsize = minimum(hgrd[shgrd[1,:]]), # Find the smallest rental unit\
        rentsupport = rentsupport, # On top of rent
        sellersupport = sellersupport, # On top of rent after selling

        ## Shared ownership
        contract_type = contract_type,
        )

    πν_array!(par,σνs,nν)

    return par
end

function SolveModel(ns)
    par=benchpar(ns=ns)
    Val,pol = solve_decproblems(par);
    aggmom, pan = simulate(par,pol,Val);

    Sol = (
        par = par,
        Val= Val,
        pol = pol,
        aggmom = aggmom,
        pan = pan,
    )
end

function πν_array!(par,σνs,nν)
    for ie in eachindex(par.egrd)
        if nν>1
            π_tmp, grd_tmp =  rouwenhorst(nν, par.ρ, σνs[ie])
            for (ia,age) in enumerate(par.agegrd)
                if age < par.Tr
                    par.πν[ia,:,:,ie] .= π_tmp
                else
                    par.πν[ia,:,:,ie] .= 0.0
                    for iν in 1:nν
                        par.πν[ia,iν,iν,ie] = 1.0
                    end
                end
            end
            par.νgrd .= exp.(grd_tmp)

            par.CDFν_ini[:,ie] = CDF_1d((par.πν[1,:,:,ie]^100)[1,:]) # Create CDF for initial ν
        else
            par.νgrd .= 1.0
            par.πν .= 1.0
            par.CDFν_ini .= 1.0
        end
    end


    for (ia,age) in enumerate(par.agegrd)
        if age < par.Tr
            if par.nϵ>1
                par.ϵgrd[:,ia] = collect(range(0.8,stop=1.2,length=par.nϵ))
            else
                par.ϵgrd[:,ia] .= 1.0
            end
        else
            iϵmid = ceil(Int,length(par.ϵgrd[:,1])/2)
            par.ϵgrd[:,ia] .= -Inf64 # Set to -Inf so if it is used everything goes to hell
            par.ϵgrd[iϵmid,ia] = 1.0 # No income shock for old retires
            par.πϵ[:,ia] .= 0.0
            par.πϵ[iϵmid,ia] = 1.0 # Only feasible outcome
        end
    end

end

function finwealth(b)
    if b > 0
        finwealth = b
    else
        finwealth = 0
    end
end

"Function to rescale all parameters to change units. E.g., move from thousands of dollars to hundreds"
function rescale_parameters(par,fac)
    par.fa .*= fac
    par.pgrd .*= fac
    par.xgrd_ini .*= fac
    par = merge(par,(xgrd = range(start=par.xgrd[1]*fac,stop=par.xgrd[end]*fac,length=length(par.xgrd)),))
    par = merge(par,(rentsupport = par.rentsupport*fac,))
    par = merge(par,(rentsupport = par.sellersupport*fac,))

    return par
end

## Real convenience function (e.g., convert to CDFs, quadratures, numerical "tricks",...)
function CDF_1d(PDFvec::Vector)
    # Check that PDFs are actually PDFs
    @assert abs(sum(PDFvec) - 1) < 0.001
    @assert minimum(PDFvec) >= 0.0
    @assert maximum(PDFvec) <= 1.0

    CDF = cumsum(PDFvec) # Running sum to get CDF

    CDF = CDF/CDF[end]
    return CDF
end

function findoutcomeCDF(pr::AbstractFloat,CDF)
    @assert pr >= 0.0 && pr <= 1.0
    @assert CDF[1] >= 0.0
    @assert CDF[end] == 1.0

    outcome = 0
    for (i,cdf) in enumerate(CDF)
        if pr <= cdf
            outcome = i
            break
        end
    end

    @assert outcome != 0
    return outcome
end

"""
    equiprob_discr(::UnivariateDistribution,::Int)

Function to discretize any UniverateNormal Distribution.
- First element is a `UnivariateDistribution` to be discretized
- Second element is an `Int` that specifies number of nodes
`grd,PDF,CDF = equiprob_discr(Normal(),11)`
"""
function equiprob_discr(dist::UnivariateDistribution,npts::Int)
    mdpnt = 1/npts*0.5
    _grd = range(mdpnt,stop=1.0 - mdpnt,length=npts)
    PMF = fill(1.0/npts,npts)
    CDF = cumsum(PMF)
    CDF ./= CDF[end]
    grd = quantile.(dist, _grd)

    return grd,PMF,CDF
end

# See https://github.com/joshday/SearchSortedNearest.jl
# Requires MIT License
# Unfortunately not released as a package
function searchsortednearest(a, x; by=identity, lt=isless, rev=false, distance=(a,b)->abs(a-b))
    i = searchsortedfirst(a, x; by, lt, rev)
    if i == 1
    elseif i > length(a)
        i = length(a)
    elseif a[i] == x
    else
        i = lt(distance(by(a[i]), by(x)), distance(by(a[i - 1]), by(x))) ? i : i - 1
    end
    return i
end

## Random/practical Functions
"Simple function to convert (s,h) to a sensible string"
function sh2string(s::AbstractFloat,h::AbstractFloat)
    @assert s ≥ 0
    @assert s ≤ 1

    if s == 0.0
        "rent (s=$(round(s;digits=1)),h="*string(h)*")"
    elseif s == 1.0
        "own (s=$(round(s;digits=1)),h="*string(h)*")"
    elseif s > 0 && s<1
        "shared (s=$(round(s;digits=1)),h="*string(h)*")"
    end
end

## Rouwenhorst algorithm, slightly modified from QuantEcon.jl. See disclaimer and copyright:
# Copyright © 2013-2019 Spencer Lyon, Thomas J. Sargent, and John Stachurski: BSD-3 All rights reserved.
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#     Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#     Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#     Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function rouwenhorst(N::Integer, ρ::Real, σ::Real, μ::Real=0.0)
    σ_y = σ / sqrt(1-ρ^2)
    p  = (1+ρ)/2
    Θ = [p 1-p; 1-p p]
    ψ = sqrt(N-1) * σ_y
    m = μ / (1 - ρ)

    state_values, p = _rouwenhorst(p, p, m, ψ, N)
    return p, state_values # Return tuple instead of QuantEcon.jl's MarkovChain object
end

function _rouwenhorst(p::Real, q::Real, m::Real, Δ::Real, n::Integer)
    if n == 2
        return [m-Δ, m+Δ],  [p 1-p; 1-q q]
    else
        _, θ_nm1 = _rouwenhorst(p, q, m, Δ, n-1)
        θN = p    *[θ_nm1 zeros(n-1, 1); zeros(1, n)] +
             (1-p)*[zeros(n-1, 1) θ_nm1; zeros(1, n)] +
             q    *[zeros(1, n); zeros(n-1, 1) θ_nm1] +
             (1-q)*[zeros(1, n); θ_nm1 zeros(n-1, 1)]

        θN[2:end-1, :] ./= 2

        return range(m-Δ, stop=m+Δ, length=n), θN
    end
end
