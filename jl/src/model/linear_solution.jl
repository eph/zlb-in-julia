include("polydef.jl")

module linear_solution

    import polydef: linsoldetails, polydetails, solutiondetails
    export get_aimsolution, lindecrule_markov
    using DSGE

function get_aimsolution(params,linsol)

    #Input
    params :: Array{Float64}
    linsol :: linsoldetails

    # Initilize Variables
    TT=Array{Float64}(42, 42)
    RR=Array{Float64}(42, 8)
    CC=Array{Float64}(42,1) 
    eu=Array{Int64}(2,1) 
    QQ=zeros(8, 8)
    GAM0=zeros(42, 42)
    GAM1=zeros(42, 42) 
    C=zeros(42,1) 
    PSI=zeros(42, 8) 
    PPI=zeros(42, 8) 

    #get aim parameters
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
    gamxhp = params[21] 
    gamdy = params[22] 
    shrgy = params[23]
    sdevtech = params[24]
    rhog = params[25] 
    sdevg = params[26] 
    rhoinv = params[27]
    sdevinv = params[28]
    rholiq = params[29] 
    sdevliq = params[30] 
    rhoint = params[31]
    sdevint = params[32]
    rhoa = params[33] 
    sdeva = params[34] 
    
    #WHAT IS GOING ON WITH THESE
    rhoelast=0
    rhoelastw=0
    sdevelast=1e-5
    sdevelastw=1e-5

    gg=1.0/(1.0-shrgy)
    gamtil=gamma/gz
    mc=(ep-1.0)/ep
    k2yrat=((mc*alpha)/(gz/beta-(1.0-delta)))*gz
    shriy=(1.0-(1.0-delta)/gz)*k2yrat
    shrcy=1.0-shrgy-shriy
    labss=( ((epw-1.0)/epw)*(1.0-alpha)*(1.0-beta*gamtil)*((ep-1.0)/ep)*(1.0/(psil*(1.0-gamtil)))*(1.0/shrcy) )^(1.0/(sigmal+1.0))
    kappaw=((1.0-gamtil)/(1.0-beta*gamtil))*epw*psil*labss^(1.0+sigmal)/phiw
    kappap=(ep-1.0)/(phi*(1.0+beta*(1.0-ap)))
    kss=labss*(gz^(alpha/(alpha-1.0)))*k2yrat^(1.0/(1.0-alpha))
    gdpss=(kss/gz)^alpha*labss^(1.0-alpha)
    invss=shriy*gdpss
    phii_jpt=phii/invss
    css=shrcy*gdpss
    rwss=(1.0-alpha)*mc*gdpss/labss
    mucss=(1.0/css)*(1.0/(1.0-gamtil))
    lamss=mucss*(1.0-beta*gamtil)
    rss=gz*pibar/beta
    rkss=gz/beta-1.0+delta

    GAM0[3, 1] = 1.0
    GAM0[2, 2] = -1.0*beta*gamtil^2.0 - 1.0
    GAM0[7, 2] = gg*shrcy
    GAM0[16, 2] = -1.0
    GAM0[21, 2] = gamtil*1.0/(-1.0*gamtil + 1.0)
    GAM0[37, 2] = -1.0
    GAM0[3, 3] = 1.0/gz*(-1.0*delta + 1.0) - 1.0
    GAM0[7, 3] = gg*shriy
    GAM0[11, 3] = phii_jpt*(beta + 1.0)
    GAM0[17, 3] = 1.0
    GAM0[22, 3] = -1.0
    GAM0[35, 3] = -1.0
    GAM0[4, 4] = -1.0
    GAM0[14, 4] = 1.0
    GAM0[15, 4] = 1.0
    GAM0[19, 4] = -1.0*kappaw
    GAM0[5, 5] = -1.0
    GAM0[9, 5] = 1.0
    GAM0[10, 5] = 1.0
    GAM0[4, 6] = -1.0
    GAM0[5, 6] = gam_dp*(-1.0*gam_rs + 1.0)
    GAM0[6, 6] = -1.0
    GAM0[18, 6] = 1.0
    GAM0[40, 6] = -1.0
    GAM0[5, 7] = gamdy*(-1.0*gam_rs + 1.0)
    GAM0[7, 7] = -1.0
    GAM0[12, 7] = 1.0
    GAM0[14, 7] = -1.0
    GAM0[5, 8] = gamxhp*(-1.0*gam_rs + 1.0)
    GAM0[8, 8] = -1.0
    GAM0[9, 9] = -1.0
    GAM0[1, 10] = -1.0
    GAM0[2, 10] = -1.0*(-1.0*gamtil + 1.0)*(-1.0*beta*gamtil + 1.0)
    GAM0[10, 10] = -1.0
    GAM0[19, 10] = -1.0*kappaw
    GAM0[42, 10] = -1.0
    GAM0[1, 11] = -1.0
    GAM0[11, 11] = -1.0
    GAM0[41, 11] = -1.0
    GAM0[8, 12] = -1.0*alpha + 1.0
    GAM0[12, 12] = alpha - 1.0
    GAM0[14, 12] = 1.0
    GAM0[15, 12] = 1.0
    GAM0[19, 12] = kappaw*sigmal
    GAM0[7, 13] = alpha*1.0/ep*gg*(ep - 1.0)
    GAM0[8, 13] = alpha
    GAM0[12, 13] = -1.0*alpha
    GAM0[13, 13] = -1.0
    GAM0[15, 13] = -1.0
    GAM0[6, 14] = kappap
    GAM0[14, 14] = -1.0
    GAM0[13, 15] = 1.0/sigmaa
    GAM0[15, 15] = -1.0
    GAM0[36, 15] = -1.0
    GAM0[16, 16] = gamtil - 1.0
    GAM0[17, 17] = -1.0
    GAM0[18, 18] = -1.0
    GAM0[19, 19] = -1.0
    GAM0[20, 19] = 1.0
    GAM0[38, 19] = -1.0
    GAM0[4, 20] = 1.0
    GAM0[20, 20] = -1.0
    GAM0[21, 21] = -1.0
    GAM0[22, 22] = -1.0
    GAM0[10, 23] = 1.0
    GAM0[27, 23] = -1.0
    GAM0[3, 24] = 1.0/gz*(-1.0*delta + 1.0) - 1.0
    GAM0[11, 24] = -1.0
    GAM0[28, 24] = -1.0
    GAM0[2, 25] = -1.0*gamtil
    GAM0[3, 25] = 1.0/gz*(-1.0*delta + 1.0)
    GAM0[4, 25] = -1.0
    GAM0[5, 25] = gamdy*(-1.0*gam_rs + 1.0)
    GAM0[11, 25] = phii_jpt
    GAM0[12, 25] = alpha
    GAM0[15, 25] = 1.0
    GAM0[16, 25] = -1.0*gamtil
    GAM0[17, 25] = 1.0
    GAM0[20, 25] = -1.0*bw + 1.0
    GAM0[29, 25] = -1.0
    GAM0[39, 25] = -1.0
    GAM0[5, 26] = 1.0
    GAM0[30, 26] = -1.0
    GAM0[7, 27] = 1.0
    GAM0[31, 27] = -1.0
    GAM0[32, 28] = -1.0
    GAM0[33, 29] = -1.0
    GAM0[12, 30] = alpha - 1.0
    GAM0[34, 30] = -1.0
    GAM0[23, 31] = -1.0
    GAM0[24, 32] = -1.0
    GAM0[25, 33] = -1.0
    GAM0[26, 34] = -1.0
    GAM0[11, 35] = -1.0*beta*phii_jpt
    GAM0[22, 35] = 1.0
    GAM0[1, 36] = -1.0*beta*1.0/gz*(-1.0*delta + 1.0) + 1.0
    GAM0[2, 37] = beta*gamtil
    GAM0[21, 37] = -1.0*1.0/(-1.0*gamtil + 1.0)
    GAM0[19, 38] = beta
    GAM0[1, 39] = -1.0
    GAM0[2, 39] = beta*gamtil
    GAM0[10, 39] = -1.0
    GAM0[11, 39] = -1.0*beta*phii_jpt
    GAM0[21, 39] = -1.0*gamtil*1.0/(-1.0*gamtil + 1.0)
    GAM0[22, 39] = 1.0
    GAM0[6, 40] = beta*1.0/(beta*(-1.0*ap + 1.0) + 1.0)
    GAM0[10, 40] = -1.0
    GAM0[1, 41] = beta*1.0/gz*(-1.0*delta + 1.0)
    GAM0[1, 42] = 1.0
    GAM0[10, 42] = 1.0

    GAM1[3, 1] = 1.0/gz*(-1.0*delta + 1.0)
    GAM1[12, 1] = alpha
    GAM1[15, 1] = 1.0
    GAM1[2, 2] = -1.0*gamtil
    GAM1[16, 2] = -1.0*gamtil
    GAM1[24, 2] = -1.0
    GAM1[11, 3] = phii_jpt
    GAM1[17, 3] = 1.0
    GAM1[25, 3] = -1.0
    GAM1[4, 4] = -1.0
    GAM1[26, 4] = -1.0
    GAM1[5, 5] = -1.0*gam_rs
    GAM1[6, 6] = -1.0*(-1.0*ap + 1.0)*1.0/(beta*(-1.0*ap + 1.0) + 1.0)
    GAM1[18, 6] = -1.0*ap + 1.0
    GAM1[20, 6] = aw - 1.0
    GAM1[5, 7] = gamdy*(-1.0*gam_rs + 1.0)
    GAM1[23, 7] = -1.0
    GAM1[27, 23] = -0.85
    GAM1[28, 24] = -1.0*rhoinv
    GAM1[30, 26] = -1.0*rhoint
    GAM1[31, 27] = -1.0*rhog
    GAM1[32, 28] = -1.0*rhoelast
    GAM1[33, 29] = -1.0*rhoelastw
    GAM1[35, 35] = -1.0
    GAM1[36, 36] = -1.0
    GAM1[37, 37] = -1.0
    GAM1[38, 38] = -1.0
    GAM1[39, 39] = -1.0
    GAM1[40, 40] = -1.0
    GAM1[41, 41] = -1.0
    GAM1[42, 42] = -1.0

    PSI[27, 1] = -1.0
    PSI[28, 2] = -1.0
    PSI[29, 3] = -1.0
    PSI[30, 4] = -1.0
    PSI[31, 5] = -1.0
    PSI[32, 6] = -1.0
    PSI[33, 7] = -1.0
    PSI[34, 8] = -1.0

    PPI[35, 1] = -1.0
    PPI[36, 2] = -1.0
    PPI[37, 3] = -1.0
    PPI[38, 4] = -1.0
    PPI[39, 5] = -1.0
    PPI[40, 6] = -1.0
    PPI[41, 7] = -1.0
    PPI[42, 8] = -1.0

    #call do_gensys(TT, CC, RR, fmat, fwt, ywt, gev, eu, loose, GAM0, GAM1, C, PSI, PPI, DIV)
    TT, CC, RR, fmat, fwt, ywt, gev, eu, loose=gensys(GAM0, GAM1, C, PSI, PPI) #DO I NEED DIV?

    QQ[1, 1] = sdevliq^2.0
    QQ[2, 2] = sdevinv^2.0
    QQ[3, 3] = sdevtech^2.0
    QQ[4, 4] = sdevint^2.0
    QQ[5, 5] = sdevg^2.0
    QQ[6, 6] = sdevelast^2.0
    QQ[7, 7] = sdevelastw^2.0
    QQ[8, 8] = sdeva^2.0

    RR = RR*sqrt.(QQ)

    linsol.pp = TT[1:28,1:28]
    linsol.sigma = RR[1:28,1:6]

    #for shocks not included in nonlinear model, set innovation part of solution matrix to zero
    if (linsol.nexogshock+linsol.nexogcont < linsol.nexog)
        linsol.sigma[:,linsol.nexogshock+1:linsol.nexog-linsol.nexogcont] = 0.0
    end

    return params,linsol
    
end



function lindecrule_markov(linsol)

    #Input
    linsol :: linsoldetails
    
    #Initilize Variables
    aalin=zeros(linsol.nvars-linsol.nexog,linsol.nvars-linsol.nexog)
    bblin=zeros(linsol.nvars-linsol.nexog,linsol.nexog)

    aalin = linsol.pp[1:linsol.nvars-linsol.nexog,1:linsol.nvars-linsol.nexog]
    for i in 1:linsol.nexogshock
        bblin[:,i] = linsol.sigma[1:linsol.nvars-linsol.nexog,i]/linsol.sigma[linsol.nvars-linsol.nexog+i,i]
    end
    
    for i in 1:linsol.nexogcont

        bblin[:,linsol.nexog-i+1] = linsol.sigma[1:linsol.nvars-linsol.nexog,linsol.nexog-i+1]/linsol.sigma[linsol.nvars-i+1,linsol.nexog-i+1]
    end
    
    return aalin,bblin

end 



end
