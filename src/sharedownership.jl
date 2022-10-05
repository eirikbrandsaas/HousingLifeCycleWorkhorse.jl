# File for deleie (s ∈ (0,1))

function housing_return_shared(par,sc,hc,pn)
    sc*pn*hc*(1.0-par.δ)
end


function ac_shared(par,st,hc,sc) # Adjustment cost (set to 0 for simplicity)
    ac = 0.0

    if st.s == sc # No change in shared ownership
        if st.h != hc # Change in size: pay purchase and sales cost
            ac += st.s*par.ms*st.p*st.h
            ac += sc*par.mb*st.p*hc
        else
            # Then there was no change, so redundant
        end
    end

    if (st.s != sc) # Change in status
        if hc == st.h # Same size, buying or selling down:
            if st.s<sc # Buying up
                ac += (sc-st.s)*par.mb*st.p*hc
            elseif st.s>sc # Selling down
                ac += (st.s - sc)*par.ms*st.p*st.h
            end
        else # Also different size, pay sales cost on what you owned plus purchase cost on new share
            ac += st.s*par.ms*st.p*st.h
            ac += sc*par.mb*st.p*hc
        end
    end

    ## Fixed (legal costs)
    ## Changing size
    # (Partial owner)      && not choose renting  && w. same size && with different s
    if (0.0 < st.s < 1.0) && (sc > 0.0) && (hc == st.h) && (sc != st.s)
        ac += par.lc
    end
    ## Selling partial ownership property
    if (0.0 < st.s < 1.0) # Current partial owner
        if sc == 0.0 # Become renter
            ac += par.ls
        elseif hc != st.h # Change size
            ac += par.ls
        end
    end

    ## Buying shared ownership
    if (0.0 < sc < 1.0) # Chose to be partial owner
        if st.s == 0.0 || st.s == 1.0 # And was either renting or owning
            ac += par.lb
        elseif hc != st.h # Or changed size
            ac += par.lb
        end
    end

    return ac
end
