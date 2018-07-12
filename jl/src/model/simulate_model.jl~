include("get_decisionrule.jl")

module simulate_model

    import polydef: polydetails, linsoldetails, solutiondetails
    import model_details: decr, decrlin
    export simulate_data, simulate_irfs


function simulate_data(capt,params,poly,linsol,alphacoeff,nonlinearswitch,seed=101293)

    # Input
    capt :: Int
    params :: Array{Float64}
    poly :: polydetails
    linsol :: linsoldetails
    alphacoeff :: Array{Float64}
    nonlinearswitch :: Bool
    seed :: Int #Optional Input, defaults to 101293

    # Initilize Variables
    msvhigh = Array{Float64}(poly.nmsv)
    msvlow = Array{Float64}(poly.nmsv)
    xrandn = Array{Float64}(poly.nexog,capt)
    modeldata=zeros(poly.nvars+2*poly.nexog,capt)
    endogvar=zeros(poly.nvars+poly.nexog,capt+1)
    innovations=zeros(poly.nexog,1)
    
    #get random normals
    iseed = MersenneTwister(seed)
    randn!(iseed,xrandn)

    counter = 0
    displacement = 400
    if (nonlinearswitch == true) 
        endogvar[1:poly.nmsv,1] = exp.(poly.endogsteady[1:poly.nmsv])
        msvhigh = exp.(poly.endogsteady[1:poly.nmsv])*2.0
        msvlow = exp.(poly.endogsteady[1:poly.nmsv])*0.01
    else
        endogvar[1:poly.nmsv,1] = poly.endogsteady[1:poly.nmsv]
        msvhigh = log.( exp.(poly.endogsteady[1:poly.nmsv])*2.0 )
        msvlow = log.( exp.(poly.endogsteady[1:poly.nmsv])*0.01 )
    end

    while true
        explosiveerror = false
        llim = counter+1
        ulim = min(capt,counter+displacement)
        for ttsim in llim:ulim
            innovations[1:poly.nexogshock] = xrandn[1:poly.nexogshock,ttsim]
            if (poly.nexogcont > 0) 
                        innovations[poly.nexog-poly.nexogcont+1:poly.nexog] = xrandn[poly.nexog-poly.nexogcont+1:poly.nexog,ttsim]
            end
            if (nonlinearswitch == true)
                endogvar[:,ttsim+1]=decr(endogvar[:,ttsim],innovations,params,poly,alphacoeff) #FIX THIS FUNCTION
                modeldata[1:poly.nvars+poly.nexog,ttsim] = endogvar[:,ttsim]
            else
                endogvar[:,ttsim+1]=decrlin(endogvar[:,ttsim],innovations,linsol) #FIX THIS FUNCTION
                modeldata[1:16,ttsim] = exp.(endogvar[1:16,ttsim+1]+poly.endogsteady[1:16])
                modeldata[20,ttsim] = exp.(endogvar[20,ttsim+1]+poly.endogsteady[20])
                modeldata[[17 18 19 21 22],ttsim] = endogvar[[17 18 19 21 22],ttsim+1]
                modeldata[poly.nvars+1:poly.nvars+poly.nexog,ttsim] = endogvar[poly.nvars+1:poly.nvars+poly.nexog,ttsim+1]
            end
            modeldata[poly.nvars+poly.nexog+1:poly.nvars+2*poly.nexog,ttsim] = innovations
        end

        for ttsim in llim:ulim
            non_explosive = ( all(msvlow.<endogvar[1:poly.nmsv,ttsim+1].<msvhigh) ) 
            explosiveerror = ( (non_explosive == false) | (isnan(endogvar[1,ttsim+1]) == true) )
            
            if (explosiveerror == true)
                counter = maximum(ttsim-350,0)
                randn!(iseed,xrandn)
                if (explosiveerror == true) 
                    println("solution exploded at ", ttsim, " vs ulim ", ulim)
                end
          
                if (counter == 0)
                    println("degenerated back to the beginning")
                end
                break
            end
        end

        if (explosiveerror == false)
            counter = counter + displacement
        end

        if (counter >= capt)
            break
        end
    end
    
    return modeldata
end 


function euler_errorsfcn(neulererrors,endogvarm1,endogvar,params,poly,alphacoeff,nlerrorswitch,linsol)

    # Input
    poly :: polydetails
    linsol :: linsoldetails
    neulererrors :: Int
    params :: Array{Float64}
    endogvarm1 :: Array{Float64}
    endogvar :: Array{Float64}
    alphacoeff :: Array{Float64}
    nlerrorswitch :: Bool

    # Initilize Variables
    euler_errorsfcn = Array{Float64}(neulererrors,1)
    endogvarp = Array{Float64}(poly.nvars+poly.nexog,1)
    lendogvar = Array{Float64}(poly.nvars+poly.nexog,1)
    lendogvarp = Array{Float64}(poly.nvars+poly.nexog,1)
    ev = Array{Float64}(6,1)
    exp_eul = zeros(6,1)
    exp_var = zeros(6,1)
    innovations = zeros(poly.nexog,1)

    beta = params[1] 
    pibar = params[2] 
    gz = params[3] 
    psil = params[4] 
    gamma = params[5] 
    sigmal = params[6] 
    phi = params[7] 
    phiw = params[8] 
    ep = params[9]  
    epw = params[10]  
    ap = params[11] 
    aw = params[12] 
    bw = params[13] 
    lamhp = params[14] 
    alpha = params[15] 
    delta = params[16] 
    phii = params[17]  
    sigmaa = params[18] 
    gam_rs = params[19] 
    gam_dp = params[20] 
    gam_xhp = params[21] 
    gam_dy = params[22] 
    shrgy = params[23] 
    rkss = params[3]/params[1]-1.0+params[16]

    for ss in 1:poly.nquad
        innovations[1:poly.nexogshock] = poly.ghnodes[:,ss] 
        if (nlerrorswitch == true)
            endogvarp = decr(endogvar,innovations,params,poly,alphacoeff)
        else
            lendogvar[1:poly.nvars] = log.(endogvar[1:poly.nvars])
            lendogvar[poly%nvars+1:poly.nvars+poly.nexog] = endogvar[poly.nvars+1:poly.nvars+poly.nexog]
            lendogvarp = decrlin(lendogvar,innovations,linsol)
            endogvarp[1:poly.nvars] = exp.(lendogvarp[1:poly.nvars])
            endogvarp[poly.nvars+1:poly.nvars+poly.nexog] = lendogvarp[poly.nvars+1:poly.nvars+poly.nexog]
        end
        techshkp = exp(endogvarp[25])
        invshkp = exp(endogvarp[24])
        ev[1] = endogvarp[10]/(endogvarp[6]*techshkp)
        utilcostp = (rkss/params[18])*(exp(params[18]*(endogvarp[13]-1.0))-1.0)
        ev[2] = (endogvarp[10]/techshkp)*( (endogvarp[15]*endogvarp[13])-utilcostp+(1.0-params[16])*endogvarp[11] )
        ev[3] = endogvarp[10]*(endogvarp[18]-1.0)*endogvarp[18]*endogvarp[7]
        ev[4] = (endogvarp[19]-1.0)*endogvarp[19] 
        ev[5] = endogvarp[16]/techshkp 
        ev[6] = endogvarp[10]*endogvarp[11]*invshkp*(endogvarp[17]-1.0)*endogvarp[17]*endogvarp[17]
        exp_var = exp_var + poly.ghweights[ss]*ev   
    end

    #ep and epw shock not incorporated, need to update if you want markup shocks
    techshk = exp(endogvar[25])
    liqshk = exp(endogvar[23])
    invshk = exp(endogvar[24])
    exp_eul[1] = (params[1]/params[3])*liqshk*endogvar[9]*exp_var[1]
    exp_eul[2] = (params[1]/params[3])*exp_var[2]/endogvar[10]
    exp_eul[3] = params[1]*exp_var[3]/(endogvar[10]*endogvar[7])+(params[9]/params[7])*( endogvar[14]-(params[9]-1.0)/params[9] ) 
    exp_eul[4] = params[1]*exp_var[4]+(params[10]/params[8])*endogvar[10]*endogvar[12]*( (params[4]*endogvar[12]^params[6]/endogvar[10])-((params[10]-1.0)/params[10])*endogvar[4] )
    exp_eul[5] = exp_var[5]
    exp_eul[6] = (params[1]/params[3])*exp_var[6]/endogvar[10] - 0.5*endogvar[11]*invshk*(endogvar[17]-1.0)*(endogvar[17]-1.0)

    euler_errorsfcn[1] = log(exp_eul[1]/endogva[10]) 
    euler_errorsfcn[2] = log(exp_eul[2]/endogvar[11])
    euler_errorsfcn[3] = exp_eul[3]-(endogvar[18]-1.0)*endogvar[18]
    euler_errorsfcn[4] = exp_eul[4]-(endogvar[19]-1.0)*endogvar[19]
    euler_errorsfcn[5] = log(exp_eul[5]/endogvar[21])
    euler_errorsfcn[6] = log( 1.0 + (1/params[18])*log(endogvar[15]/rkss) )-log(endogvar[13])
    euler_errorsfcn[7] = endogvar[1]-(1.0-delta)*(endogvarm1[1]/(gz*techshk)) -invshk*endogvar[3]*(1.0-(phii/2.0)*(endogvar[17]-1.0)^2) #capital acc
    euler_errorsfcn[8] = endogvar[2]-(gamma*endogvarm1[2]/(gz*techshk)+1.0/endogvar[16])  #defines cc using muc
    euler_errorsfcn[9] = endogvar[11]*invshk*(endogvar[17]-1.0)*endogvar[17] + (1.0/phii)*(1.0-endogvar[11]*invshk)-exp_eul[6] #investment demand

    dwtildem1 = (pibar^aw)*(endogvarm1[6]^(1.0-aw))
    gzwage = gz*techshk^(1.0-bw)
    dw = endogvar[19]*dwtildem1*gzwage
    euler_errorsfcn[10] = endogvar[4]-(endogvarm1[4]*dw)/(gz*techshk*endogvar[6])

    xhp = alpha*log(endogvar[13]) + (1-alpha)*(log(endogvar[12])-poly.endogsteady[12])
    rss = gz*pibar/beta
    euler_errorsfcn[11] = log(rss)+gam_rs*(log(endogvarm1[5])-log(rss)) + (1.0-gam_rs)*( gam_dp*(log(endogvar[6]/pibar)) +
                            gam_dy*log(endogvar[7]*techshk/endogvarm1[7]) + gam_xhp*xhp ) + endogvar[26] - log(endogvar[5])

    vp = endogvar[18]  
    dptildem1 = (pibar^ap)*(endogvarm1[6]^(1.0-ap))
    euler_errorsfcn[12] = vp*dptildem1 - endogvar[6]

    gss = 1.0/(1.0-shrgy)
    gshk = exp( log(gss) + endogvar[27] )
    aayy = 1.0/gshk-(phi/2.0)*(vp-1.0)^2
    utilcost = (rkss/sigmaa)*(exp(sigmaa*(endogvar[13]-1.0))-1.0)
    euler_errorsfcn[13] = (1.0/aayy)*( endogvar[2]+endogvar[3] + utilcost*(endogvarm1[1]/(gz*techshk)) ) - endogvar[7]

    euler_errorsfcn[14] = endogvar[7]^(1.0/(1.0-alpha))*(endogvar[13]*endogvarm1[1]/(gz*techshk))^(alpha/(alpha-1.0)) - endogvar[12] #labor
    euler_errorsfcn[15] = endogvar[14]-endogvar[4]*endogvar[12]/((1.0-alpha)*endogvar[7])  #mc
    euler_errorsfcn[16] = (alpha/(1.0-alpha))*(endogvar[4]*endogvar[12]*gz*techshk/(endogvar[13]*endogvarm1[1])) - endogvar[15] #rentalk
    euler_errorsfcn[17] = endogvar[10] + (gamma/gz)*beta*endogvar[21]  - endogvar[16] #muc

    #average absolute value of euler errors
    for i in 1:neulererrors-1
        euler_errorsfcn[i] = abs(euler_errorsfcn[i])
        euler_errorsfcn[neulererrors] = euler_errorsfcn[neulererrors] + euler_errorsfcn[i]/17
    end
    
    return euler_errorsfcn

end




function simulate_irfs(capt,ndraw,params,poly,alphacoeff,linsol,endogvarbas0,endogvarshk0,innov0,neulererrors)

    # Input
    capt :: Int
    ndraw :: Int
    neulererrors :: Int
    poly :: polydetails
    linsol :: linsoldetails
    alphacoeff :: Array{Float64}
    params :: Array{Float64}
    endogvarbas0 :: Array{Float64}
    endogvarshk0 :: Array{Float64}
    innov0 :: Array{Float64}
          
    # Initilize variables
    endogvarbas = Array{Float64}(poly.nvars+poly.nexog,capt+1)
    endogvarshk = Array{Float64}(poly.nvars+poly.nexog,capt+1)
    endoglinbas = Array{Float64}(poly.nvars+poly.nexog,capt+1)
    endoglinshk = Array{Float64}(poly.nvars+poly.nexog,capt+1)
    endoglinshkm1_exp = Array{Float64}(poly.nvars+poly.nexog,1)
    endoglinshk_exp = Array{Float64}(poly.nvars+poly.nexog,1)
    varsort = Array{Float64}(ndraw,1)
    innovplus = Array{Float64}(poly.nexog,1)
    xrandn = Array{Float64}(poly.nexog,capt*ndraw)
    errors_nl = Array{Float64}(neulererrors,1)
    errors_lin = Array{Float64}(neulererrors,1)
    wishlist_level = Array{Float64}(Int,nwish_level)   
    endogirf = zeros(poly.nvars+poly.nexog+2,capt)  #+2 should match nwish_level
    linirf = zeros(poly.nvars+poly.nexog+2,capt)  #+2 shoud match nwish_level
    euler_errors = zeros(2*neulererrors,capt)  
    levelresponse = zeros(nwish_level*capt,ndraw)
    levellin = zeros(nwish_level*capt,ndraw)
    innovations = zeros(poly.nexog,1)
    premirf = zeros(2,capt)
    premiumbas = ones(capt,1)
    premiumshk = ones(capt,1)
    nwish_level = 2
    eulererrorswitch = true    
    seed = 101293

    #get random normals
    iseed = MersenneTwister(seed)
    randn!(iseed,xrandn)

    #9=nominal rate, 5=notional rate
    wishlist_level = [9 5]
    j_used = 0
    for j in 1:ndraw
        endogvarbas[:,1] = endogvarbas0
        endogvarshk[:,1] = endogvarshk0
        endoglinbas[1:poly.nvars,1] = log.(endogvarbas0[1:poly.nvars])
        endoglinbas[poly.nvars+1:poly.nvars+poly.nexog,1] = endogvarbas0[poly.nvars+1:poly.nvars+poly.nexog]
        endoglinshk[1:poly.nvars,1] = log.(endogvarshk0[1:poly.nvars])
        endoglinshk[poly.nvars+1:poly.nvars+poly.nexog,1] = endogvarshk0[poly.nvars+1:poly.nvars+poly.nexog]
        for ttsim in 1:capt   
            innovations[1:poly.nexogshock] = xrandn[1:poly.nexogshock, capt*(j-1)+ttsim]
            if (poly.nexogcont > 0) 
                    innovations[poly.nexog-poly.nexogcont+1:poly.nexog] = xrandn[poly.nexog-poly.nexogcont+1:poly.nexog,ttsim]
            end
            endogvarbas[:,ttsim+1] = decr(endogvarbas[:,ttsim],innovations,params,poly,alphacoeff)
            endoglinbas[:,ttsim+1] = decrlin(endoglinbas[:,ttsim],innovations,linsol)
            if (ttsim == 1)
                innovplus = innovations + innov0
            else
                innovplus = innovations
            end
            endogvarshk[:,ttsim+1] = decr(endogvarshk[:,ttsim],innovplus,params,poly,alphacoeff) 
            endoglinshk[:,ttsim+1] = decrlin(endoglinshk[:,ttsim],innovplus,linsol)

            #calc_premium(endogvarshk(:,ttsim),params,poly,alphacoeff,premiumshk(ttsim))
            #calc_premium(endogvarbas(:,ttsim),params,poly,alphacoeff,premiumbas(ttsim))
      
            #IRFs (any([nquadsingle == i for i in nquadsingleset]) == false)
            if any([isnan(i) for i in endogvarshk[:,ttsim+1]]) | any([isnan(i) for i in endogvarbas[:,ttsim+1]])
                println("failure at draw = ", j)
                break
            else
         
                #compute absolute value of euler errors
                if (eulererrorswitch == true) 
                    nlerrorswitch = true
                    errors_nl = euler_errorsfcn(neulererrors,endogvarshk[:,ttsim],endogvarshk[:,ttsim+1],params,poly,alphacoeff,nlerrorswitch,linsol)
                    if any([isnan(i) for i in errors_nl])
                        break
                    else 
                        euler_errors[1:neulererrors,ttsim] = euler_errors[1:neulererrors,ttsim] + errors_nl
                    end
            
                    nlerrorswitch = false
                    endoglinshkm1_exp[1:poly.nvars] = exp.(endoglinshk[1:poly.nvars,ttsim])
                    endoglinshkm1_exp[poly.nvars+1:poly.nvars+poly.nexog] = endoglinshk[poly.nvars+1:poly.nvars+poly.nexog,ttsim]
                    endoglinshk_exp[1:poly.nvars] = exp.(endoglinshk[1:poly.nvars,ttsim+1])
                    endoglinshk_exp[poly.nvars+1:poly.nvars+poly.nexog] = endoglinshk[poly.nvars+1:poly.nvars+poly.nexog,ttsim+1]
                    errors_lin = euler_errorsfcn(neulererrors,endoglinshkm1_exp,endoglinshk_exp,params,poly,alphacoeff,nlerrorswitch,linsol)
                    if any([isnan(i) for i in errors_lin])
                        break
                    else 
                        euler_errors[neulererrors+1:2*neulererrors,ttsim] = euler_errors[neulererrors+1:2*neulererrors,ttsim] + errors_lin   
                    end
                end

                if (ttsim == capt)
                    j_used = j_used + 1
                end
         
                endogirf[1:poly.nvars+poly.nexog,ttsim+1] = endogirf[1:poly.nvars+poly.nexog,ttsim+1] + 100.0*log.(endogvarshk[:,ttsim+1]/endogvarbas[:,ttsim+1])  
                linirf[1:poly.nvars+poly.nexog,ttsim] = linirf[1:poly.nvars+poly.nexog,ttsim] + 100.0*(endoglinshk[:,ttsim+1]-endoglinbas[:,ttsim+1]) 
                premirf[1,ttsim] = premirf[1,ttsim] + 100.0*log(premiumshk[ttsim]/premiumbas[ttsim])

                for i in 1:nwish_level   
                    levelresponse[(i-1)*capt+ttsim,j_used] = 400.0*log(endogvarshk[wishlist_level[i],ttsim+1])
                    levellin[(i-1)*capt+ttsim,j_used] = 400.0*endoglinshk[wishlist_level[i],ttsim+1]
                end
            end
        end
    end

    endogirf = endogirf/j_used
    linirf = linirf/j_used
    premirf = premirf/j_used
    euler_errors = euler_errors/j_used

    return endogirf,linirf,euler_errors,premirf
                                                        
end 
                            
end
