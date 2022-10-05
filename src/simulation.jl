function initial_state(pr,par)
    xs = par.xgrd_ini[findoutcomeCDF(pr.x_ini,par.CDFx_ini)]
    ie = findoutcomeCDF(pr.e_ini,par.CDFe_ini)
    iν = findoutcomeCDF(pr.ν_ini,view(par.CDFν_ini,:,ie))
    is = 1
    ih = 1
    ip = ceil(Int,length(par.pgrd)/2)
    ia = 1

    return xs,ie,iν,is,ih,ip,ia
end

function indshocks(par)
    rng = StableRNG(par.seed)
    len = length(par.agegrd)
    indprob = [(
            ν = rand(rng,len),
            p = rand(rng,len),
            ν_ini = rand(rng,),
            x_ini = rand(rng,),
            e_ini = rand(rng,),
    ) for i=1:par.N]


    return indprob
end

function simpar(par) # Create numerical objects for simulation
    CDFν = similar(par.πν);
    for ie in eachindex(par.egrd)
        for ia in eachindex(par.agegrd)
            for iν in 1:length(par.νgrd)
                CDFν[ia,iν,:,ie] = CDF_1d(par.πν[ia,iν,:,ie])
            end
        end
    end

    CDFp = similar(par.πp)
    for ip in eachindex(par.pgrd)
        CDFp[ip,:] = CDF_1d(par.πp[ip,:])
    end


    spar = (CDFν = CDFν,
            CDFp = CDFp)
end


function simulate_ind!(par,spar,probs,Vcond,ccond,cvec,bvec,svec,hvec,xvec,pvec,νvec,evec,fvec)
    xs, ie, iν, is, ih, ip, ia = initial_state(probs,par)
    νvec[ia] = iν
    svec[ia] = is
    hvec[ia] = ih
    pvec[ia] = ip
    xvec[ia] = xs
    evec .= ie

    _tmp = fill(-Inf64,(length(par.sgrd),length(par.hgrd)))

    for (ia,age) in enumerate(par.agegrd)
        # Need to be careful with interpolation of discrete choices
        _tmp .= -Inf64
        for ihn in eachindex(par.hgrd)
            for isn in eachindex(par.sgrd)
                if par.shgrd[isn,ihn] == true # Problem only exists for possible combination of s,h
                    _tmp[isn,ihn] = Vcond[isn,ihn,ie,iν,is,ih,ip,ia](xs)
                end
            end
        end
        (isc,ihc) = Tuple(argmax(_tmp))

        sc = par.sgrd[isc]
        hc = par.hgrd[ihc]
        cc = ccond[isc,ihc,ie,iν,is,ih,ip,ia](xs)
        st = state(par,xs,ie,par.hgrd[ih],par.sgrd[is],iν,ip,par.pgrd[ip],ia)
        bc = bchoice(par,st, cc,hc,sc)
        cvec[ia] = cc
        bvec[ia] = bc
        fvec[ia] = feasible(par,st,cc,bc,sc,hc)

        if age < par.Te
            νvec[ia+1] = iν = findoutcomeCDF(probs.ν[ia+1],view(spar.CDFν,ia+1,iν,:,ie))
            svec[ia+1] = is = isc
            hvec[ia+1] = ih = ihc
            pvec[ia+1] = ip = findoutcomeCDF(probs.p[ia+1],view(spar.CDFp,ip,:))
            xvec[ia+1] = xs = LoM(par,bc,sc,hc,par.pgrd[ip])
        end
    end
end

function simulate(par,pol,Val)
    # Set up necessary objects
    spar = simpar(par)
    probs = indshocks(par) # Draw shocks for all households
    _Vtmp = replace!(copy(Val.Vcond),-Inf=>nextfloat(typemin(Val.Vcond[1]))) # Interpolation algorithm get's messed up by the -Inf, so replace with the most negative non-inf Float
    Vcond_itp = [linear_interpolation(par.xgrd,_Vtmp[isc,ihc,:,ie,iν,is,ih,ip,ia],extrapolation_bc=Line()) for isc in eachindex(par.sgrd),ihc in eachindex(par.hgrd),  ie=1:length(par.egrd),iν=1:length(par.νgrd),is=1:length(par.sgrd),ih=1:length(par.hgrd),ip=1:length(par.pgrd),ia in eachindex(par.agegrd)]
    Ccond_itp = [linear_interpolation(par.xgrd,pol.ccond[isc,ihc,:,ie,iν,is,ih,ip,ia],extrapolation_bc=Line()) for isc in eachindex(par.sgrd),ihc in eachindex(par.hgrd),  ie=1:length(par.egrd),iν=1:length(par.νgrd),is=1:length(par.sgrd),ih=1:length(par.hgrd),ip=1:length(par.pgrd),ia in eachindex(par.agegrd)]

    # Preallocate matrices for simulation (each column is an individual)
    cmat = fill(-Inf64,length(par.agegrd),par.N)
    bmat = fill(-Inf64,length(par.agegrd),par.N)
    xmat = fill(-Inf64,length(par.agegrd),par.N)
    ismat = fill(0,length(par.agegrd),par.N)
    ihmat = fill(0,length(par.agegrd),par.N)
    ipmat = fill(0,length(par.agegrd),par.N)
    iνmat = fill(0,length(par.agegrd),par.N)
    iemat = fill(0,length(par.agegrd),par.N)
    ifmat = fill(true,length(par.agegrd),par.N) # = true if choice is feasible

    # Simulate all households
    @Threads.threads for id = 1:par.N
        @views simulate_ind!(par,spar,probs[id],Vcond_itp,Ccond_itp,cmat[:,id],bmat[:,id],ismat[:,id],ihmat[:,id],xmat[:,id],ipmat[:,id],iνmat[:,id],iemat[:,id],ifmat[:,id])
    end

    # Calculate moments by age
    agg_mom = (
        age = par.agegrd,
        x = mean(xmat;dims=2),xse = std(xmat;dims=2),
        c = mean(cmat;dims=2),cse = std(cmat;dims=2),
        b = mean(bmat;dims=2),bse = std(bmat;dims=2),
        rent = mean(par.sgrd[ismat].==0.0;dims=2),
        own = mean(par.sgrd[ismat].==1.0;dims=2),
        h = mean(par.hgrd[ihmat];dims=2),
        iν = mean(iνmat;dims=2),
        ip = mean(ipmat;dims=2),
        ie = mean(iemat;dims=2),
        feasible = mean(ifmat;dims=2),
    )

    # Create a Panel (stored in a NamedTuple)
    pan = (
        id = sort!(repeat(1:1:par.N,length(par.agegrd))),
        age = repeat(par.agegrd,par.N),
        c = vec(cmat),
        b = vec(bmat),
        x = vec(xmat),
        is = vec(ismat),
        ih = vec(ihmat),
        ip = vec(ipmat),
        iν = vec(iνmat),
        ie = vec(iemat),
        feasible = vec(ifmat),
        )

    return agg_mom, pan
end
