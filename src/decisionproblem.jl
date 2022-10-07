function util_exp(par,st,cc,hc,sc,Vnxt)
    bc = bchoice(par,st,cc,hc,sc)
    val = util(cc,hc,sc,par)
    for iϵn in 1:par.nϵ
        for ipn in eachindex(par.pgrd)
            for iνn in eachindex(par.νgrd)
                prob = par.πν[st.ia,st.iν,iνn,st.ie]*par.πp[st.ip,ipn]*par.πϵ[iϵn,st.ia]
                if prob > 0 # Only calculate next-period val if prob>0
                    xn = LoM(par,bc,sc,hc,par.pgrd[ipn]) + trans_inc_shock(par,st,iνn,par.ϵgrd[iϵn,st.ia])
                    val += par.β*prob*Vnxt[iνn,ipn](xn)
                end
            end
        end
    end

    return val
end

function decproblem(par,st,Vnxt,vcond,ccond)
    for (ih,hc) in enumerate(par.hgrd)
        for (is,sc) in enumerate(par.sgrd)
            if feasible_discrete(par,is,ih,st,sc) == true  # Some combination of discrete choices can be ruled out by assumption
                lb = minc()
                ub = maxc(par,st,hc,sc)
                if ub ≥ lb # Otherwise no feasible choice
                    res = optimize(x -> -util_exp(par,st,x,hc, sc,view(Vnxt,:,:,is,ih,st.ie)), lb, ub)
                    vcond[is,ih] = -res.minimum
                    ccond[is,ih] = res.minimizer
                end
            end
        end
  end
  pol = argmax(vcond)
  is = pol[1]
  ih = pol[2]

  return vcond[pol], is, ih, ccond[pol]
end

function solve_decproblems(par)
  stspace = length(par.xgrd),length(par.egrd),length(par.νgrd),length(par.sgrd),length(par.hgrd),length(par.pgrd),length(par.agegrd)
  ## Initialize elements
  ihc = fill(0,stspace)
  isc = fill(0,stspace)
  cc = fill(-Inf64,stspace)
  bc = fill(-Inf64,stspace)
  V = fill(-Inf64,stspace .+ (0,0,0,0,0,0,1)) # (Note one extra element for age)
  Vcond = fill(-Inf64,(length(par.sgrd),length(par.hgrd),stspace...))
  Ccond = fill(-Inf64,(length(par.sgrd),length(par.hgrd),stspace...))

  V[:,:,:,:,:,:,end]  .= 0.0 # Set value to 0 after death

  for ia in reverse(eachindex(par.agegrd))
    Vnxt = [linear_interpolation(par.xgrd,V[:,ie,iν,is,ih,ip,ia+1],extrapolation_bc=Line())
        for iν=1:length(par.νgrd),ip=1:length(par.pgrd),is=1:length(par.sgrd),ih=1:length(par.hgrd),ie=1:length(par.egrd)]
    @Threads.threads for (ix,xs) in collect(enumerate(par.xgrd))
      for (ie,es) in enumerate(par.egrd)
        for (iν,νs) in enumerate(par.νgrd)
          for (is,ss) in enumerate(par.sgrd)
            for (ih,hs) in enumerate(par.hgrd)
              for (ip,ps) in enumerate(par.pgrd)
                  if par.shgrd[is,ih] == true # Problem only exists for possible combination of s,h
                      st = state(par,xs,ie,hs,ss,iν,ip,ps,ia)
                      ist = LinearIndices(stspace)[ix,ie,iν,is,ih,ip,ia]
                      V[ist], isc[ist], ihc[ist],cc[ist] = decproblem(par,st,Vnxt,view(Vcond,:,:,ix,ie,iν,is,ih,ip,ia),view(Ccond,:,:,ix,ie,iν,is,ih,ip,ia))
                      bc[ist] = bchoice(par,st, cc[ist],par.hgrd[ihc[ist]],par.sgrd[isc[ist]])

                      #Throws assertion error if choice is feasible
                      assert_feasible(par,st,cc[ist],bc[ist],par.sgrd[isc[ist]],par.hgrd[ihc[ist]])
                  end
              end
            end
          end
        end
      end
    end
  end

  Val = (V=V,Vcond=Vcond)
  pol = (c=cc,b=bc,is=isc,ih=ihc,ccond=Ccond)
  return Val, pol
end
