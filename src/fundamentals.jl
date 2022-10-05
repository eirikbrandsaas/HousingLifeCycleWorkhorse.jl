function state(par,x,ie,h,s,iν,ip,p,ia)
    @assert s ≤ 1.0
    @assert h ≥ 1.0
    st = (x = x,
        ie = ie,
        h = h,
        s = s,
        iν = iν,
        y = wage(par,ia,iν,ie,p,x,h,s),
        ip = ip,
        p = p,
        ia = ia)
end

function util(c,h,s,par)
    ((c^(1.0-par.η)*h^(par.η)*ownshift(s,par))^(1.0-par.γ)-1.0)/(1.0-par.γ)
end

function ownshift(s,par)
    1.0+(s^par.α)*par.χ
end

function bchoice(par,st,cc,hc,sc)
    b = st.x + st.y - cc - ac(par,st,hc,sc) - housing_exp(par,st.p,hc,sc)
end

function ac(par,st,hc,sc) # Adjustment cost (set to 0 for simplicity)
    ac = 0.0
    if (0.0 < st.s < 1.0) || (0.0 < sc < 1.0) # If housing status state or choice is ∉{0,1}, use shared ownership function
        ac += ac_shared(par,st,hc,sc)
    else # Standard model with only sales/purchases cost
        if (st.s == 1.0) && (st.h != hc || sc == 0.0)
            ac += par.ms*st.p*st.h
        end
        if (sc == 1.0) && (st.h != hc || st.s == 0.0)
            ac += par.mb*st.p*hc
        end
    end
    return ac
end

function housing_exp(par,ps,hc,sc)
    ps*((1-sc)*par.κ*hc + sc*hc)
end

function wage(par,ia,iν,ie,ps,xs,hs,ss)
    w = par.νgrd[iν]*par.fa[ia,ie]
    w = max(w,minwage(par,ps))
    w = sellersupport(par,w,xs,hs,ss,ps)

    return w
end

function minwage(par,ps)
    par.minrentalsize*ps*par.κ + par.rentsupport
end

function sellersupport(par,ws,xs,hs,ss,ps)
    if ss == 1.0
        mincost_selling_to_rent = par.ms*ps*hs*ss + par.κ*ps*par.minrentalsize
        if ws  < mincost_selling_to_rent - xs
            ws = mincost_selling_to_rent - xs + par.sellersupport
        end
    end

    return ws
end
function maxborr(par,st,hc::Float64,sc::Float64)
    min(BTI(par,sc,st.y),LTV(par,hc,sc,st.p))
end

function LTV(par,hc,sc,ps)
    (1.0-par.d)*ps*hc*sc
end

function BTI(par,sc,ws)
    par.b2y*sc*ws
end

function minc() # Minimum consumption
    0.0
end

function maxc(par,st,hc,sc) # Maximum consumption
    st.x + st.y - ac(par,st,hc,sc) + maxborr(par,st,hc,sc) - housing_exp(par,st.p,hc,sc)
end

function assert_feasible(par,st,cc,bc,sc,hc)
    @assert cc ≥ minc() "$st, \n c=$cc, b=$bc ≤ minc = $(minc())"
    @assert cc ≤ maxc(par,st,hc,sc)  "$st, c=$cc, maxc = $(maxc(par,st,hc,sc))"
    @assert bc ≥ -maxborr(par,st,hc,sc) "$st, c=$cc,b=$bc ≤ $(-maxborr(par,st,hc,sc))"
end

function feasible(par,st,cc,bc,sc,hc)
    f = true
    if cc < minc()
            f=false
    end
    if cc > maxc(par,st,hc,sc)
        f = false
    end
    if bc< -maxborr(par,st,hc,sc)
        f = false
    end

    return f
end

function LoM(par,bc,sc,hc,pn)
    xn = bc*(1.0+r(par,bc))
    if sc == 0.0
        # nothing
    elseif sc == 1.0
        xn += pn*hc*(1.0-par.δ)
    else
        xn += housing_return_shared(par,sc,hc,pn)
    end
    return xn
end

function r(par,bc)
    par.rf + par.rm*(bc<0.0)
end

function feasible_discrete(par,is,ih,st,sc)
    feasible = false
    if par.shgrd[is,ih] == true
        feasible = true
    end

    return feasible
end
