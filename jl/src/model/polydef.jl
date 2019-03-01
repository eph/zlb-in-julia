module polydef
    import FastGaussQuadrature: gausshermite
    export linsoldetails, polydetails, solutiondetails, initializesolution!, setgridsize, 
        exoggridindex, ghquadrature,sparsegrid, smolyakpoly, initializelinearsolution!, initializetestsolution!, gensys
using LinearAlgebra

mutable struct linsoldetails
    
    nparams :: Int 
    nvars :: Int
    nexog :: Int 
    nexogshock :: Int
    nexogcont :: Int
    pp :: Array{Float64} 
    sigma :: Array{Float64} 
    endogsteady :: Array{Float64}
    linsoldetails()=new() #Initilizes all the values to undefined values, closet approximation to allocatable array
    
end 

mutable struct polydetails
    
    nfunc :: Int
    nmsv :: Int
    nvars :: Int
    ngrid :: Int
    nparams :: Int
    nindplus :: Int
    nexog :: Int
    nexogshock :: Int
    nexogcont :: Int
    ns :: Int
    nsexog :: Int
    ninter :: Int
    nquad :: Int
    zlbswitch :: Bool
    indplus :: Array{Int}
    nshockgrid :: Array{Int}
    interpolatemat :: Array{Int}
    endogsteady :: Array{Float64}
    slopeconmsv :: Array{Float64}
    shockbounds :: Array{Float64}
    shockdistance :: Array{Float64}
    exoggrid :: Array{Float64}
    ghnodes :: Array{Float64}
    ghweights :: Array{Float64}
    polydetails()=new() #Initilizes all the values to undefined values, closet approximation to allocatable array
    
end

mutable struct solutiondetails
    
    poly :: polydetails
    linsol :: linsoldetails
    number_shock_values :: Int64
    startingguess :: Bool
    exogvarinfo :: Array{Int64}
    bbt :: Array{Float64}
    bbtinv :: Array{Float64}
    xgrid :: Array{Float64}
    alphacoeff :: Array{Float64}
    solutiondetails()=new()
    
end

function setgridsize(nexog,nshockgrid)

    #Input
    nexog :: Int
    nshockgrid :: Array{Int}

    nexogshock = 0
    for i in 1:nexog
        if (nshockgrid[i] > 1)
            nexogshock = nexogshock + 1
        else
            break
        end
    end
    ninter = 2^nexogshock

    #allocate matrices for shock processes
    ns = 1
    number_shock_values = 0
    for i in 1:nexogshock
        ns = ns*nshockgrid[i]
        number_shock_values = number_shock_values + nshockgrid[i]   
    end
    
    return nexogshock,ninter,ns,number_shock_values

end 

function exoggridindex(ngrid,nexog,ns)
      
    # Input
    nexog :: Int
    ns :: Int
    ngrid :: Array{Int}

    #Initilize Variables
    exoggridindex = zeros(Int,nexog,ns)
    blocksize = 1 #Initilize blocksize
    
    for ie in nexog:-1:1
        if ie == nexog
            blocksize = 1
        else
            blocksize = blocksize*ngrid[ie+1]
        end
        ncall = div(ns,blocksize*ngrid[ie]) #use div function to keep Int type
        for ic in 1:ncall
            for ib in 1:ngrid[ie]
                exoggridindex[ie,ngrid[ie]*blocksize*(ic-1)+blocksize*(ib-1)+1:ngrid[ie]*blocksize*(ic-1)+blocksize*ib] .= ib
            end
        end
    end
    
    return exoggridindex
    
end 

function ghquadrature(nquadsingle,nexog) 

    # Input
    nquadsingle :: Int64
    nexog :: Int64

    # Initilize Variables
    quadnodes_s=zeros(nquadsingle,1)
    quadweights_s=zeros(nquadsingle,1)
    ghnodes=zeros(nexog,nquadsingle^nexog)#Not sure I should make this zeros
    ghweights_mat=zeros(nexog,nquadsingle^nexog)#Not sure I should make this zeros
    ghweights=Array{Float64}(undef,nquadsingle^nexog,1)#Not sure I should make this zeros
    #const const_pi = 3.14159265358979323846
    const_pi = 3.14159265358979323846
    
    quadnodes_s,quadweights_s=gausshermite(nquadsingle) ##

    nquad = nquadsingle^nexog
    blocksize = 1 #Must Initilize blocksize
    for ie = nexog:-1:1
        if (ie == nexog)
            blocksize = 1
        else
            blocksize = blocksize*nquadsingle # it was nquadsingle*blocksize ??
        end
        ncall = div(nquad,nquadsingle*blocksize)
        for ic in 1:ncall
            for ib in 1:nquadsingle
                left=nquadsingle*blocksize*(ic-1)+blocksize*(ib-1)+1
                right=nquadsingle*blocksize*(ic-1)+blocksize*ib
                ghnodes[ie,left:right] .= sqrt(2)*quadnodes_s[ib]
                ghweights_mat[ie,left:right] .= quadweights_s[ib]
            end
        end
    end

    ghweights = (1.0/const_pi)^(nexog/2.0)*prod(ghweights_mat,dims=1) # product of ghweights_mat along the first dimension

    return nquad,ghnodes,ghweights
    
end 

function smolyakpoly(nmsv,ngrid,nindplus,indplus,xx)

    # Input
    nmsv :: Int
    ngrid :: Int
    nindplus :: Int
    indplus :: Array{Int}
    xx:: Array{Float64}
    
    # Initilize Variables
    smolyakpoly=Array{Float64}(undef,ngrid,1)

    smolyakpoly[1] = 1.0
    for i in 1:nmsv
	smolyakpoly_aux = xx[i]
        smolyakpoly[2*i] = smolyakpoly_aux 
        smolyakpoly[2*i+1] = 2.0*(smolyakpoly_aux)^2-1.0
    end

    for i in 1:nindplus
	xx_aux =xx[indplus[i]]
        smolyakpoly[2*nmsv+2*(i-1)+2] = 4.0*xx_aux^3-3.0*xx_aux
        smolyakpoly[2*nmsv+2*(i-1)+3] = 8.0*xx_aux^4-8.0*xx_aux^2+1.0
    end
    
    return smolyakpoly
    
end 


function sparsegrid(nmsv,nindplus,ngrid,indplus)

    # Input 
    nmsv :: Int
    nindplus :: Int
    ngrid :: Int
    indplus :: Array{Int}

    #Initilize Variables
    xgrid = zeros(nmsv,ngrid)
    bbt = zeros(ngrid,ngrid)
    
    for i in 1:nmsv
        xgrid[i,2*i] = -1.0
        xgrid[i,2*i+1] = 1.0
    end

    for i in 1:nindplus
        xgrid[indplus[i],2*nmsv+2*(i-1)+2] = -1.0/sqrt(2.0)
        xgrid[indplus[i],2*nmsv+2*(i-1)+3] = 1.0/sqrt(2.0)
    end

    #form bbt matrix 
    for i in 1:ngrid
	bbt_aux = smolyakpoly(nmsv,ngrid,nindplus,indplus,xgrid[:,i])
        bbt[:,i] = bbt_aux
    end

    
    # find bbt inverse matrix
    bbtinv = copy(bbt)
    bbtinv,ipiv,info=LinearAlgebra.LAPACK.getrf!(bbtinv)
    if (info == 0) 
        LinearAlgebra.LAPACK.getri!(bbtinv,ipiv)
    else
        println("something went wrong with getrf! (sparsegrid)")
        println("info = ", info)
    end
    
    return xgrid,bbt,bbtinv

end

function initializelinearsolution!(nparams,nvars,nexog,nexogshock,nexogcont,linsol)  

    # Input
    nparams :: Int
    nvars :: Int
    nexog :: Int
    nexogshock :: Int
    nexogcont :: Int
    linsol :: linsoldetails

    linsol.nparams = nparams
    linsol.nvars = nvars+nexog  #nvars is equal to poly%nvars which excludes exogenous shocks
    linsol.nexog = nexog
    linsol.nexogshock = nexogshock
    linsol.nexogcont = nexogcont
    
    linsol.pp=Array{Float64}(undef,linsol.nvars,linsol.nvars)
    linsol.sigma=Array{Float64}(undef,linsol.nvars,linsol.nexog)
    linsol.endogsteady=Array{Float64}(undef,linsol.nvars,1)
        
    return linsol

end 

function initializesolution!(solution) # ! to indicate that this function mutates and input

    #I don't think we need to declare types of these in the future (may not even work)
    nquadsingle =3
    
    #put shocks into polynomial approximation if necessary
    nexogadj = solution.poly.nexog - solution.poly.nexogcont
    nmsvadj = solution.poly.nmsv + solution.poly.nexogcont

    solution.poly.ngrid = 2*(solution.poly.nmsv+solution.poly.nexogcont)+2*solution.poly.nindplus+1

    #set nexogshock,ns, and number_shock_values
    solution.poly.nexogshock,solution.poly.ninter,solution.poly.ns,solution.number_shock_values=setgridsize(nexogadj,solution.poly.nshockgrid)
    solution.poly.nquad = nquadsingle^(solution.poly.nexogshock+solution.poly.nexogcont)

    #Set exogvarinfo
    solution.exogvarinfo=Array{Int64}(undef,nexogadj,solution.poly.ns)
    solution.exogvarinfo[1:solution.poly.nexogshock,:] = exoggridindex(solution.poly.nshockgrid,solution.poly.nexogshock,solution.poly.ns)
    solution.exogvarinfo[solution.poly.nexogshock+1:nexogadj,:] .= 1

    #get matrix used for interpolating the shocks
    solution.poly.interpolatemat=Array{Int64}(undef,solution.poly.nexogshock,2^solution.poly.nexogshock)
    blocksize = 1
    for i in solution.poly.nexogshock:-1:1
        if (i == solution.poly.nexogshock)
            blocksize = 1
        else
            blocksize = 2*blocksize
        end
        #blocksize=2^(solution.poly.nexogshock-i)
        ncall = div(2^(solution.poly.nexogshock-1),blocksize)
        for j in 1:ncall
            for k = 1:2
                left=2*blocksize*(j-1)+blocksize*(k-1)+1
                right=2*blocksize*(j-1)+blocksize*k
                solution.poly.interpolatemat[i,left:right] .= k-1 
            end
        end
    end

    #get quadrature nodes and weights
    nquadadj = solution.poly.nexogshock+solution.poly.nexogcont
    solution.poly.nquad,solution.poly.ghnodes,solution.poly.ghweights=ghquadrature(nquadsingle,nquadadj)

    #construct sparse grid, bb matrix and its inverse
    solution.xgrid,solution.bbt,solution.bbtinv=sparsegrid(nmsvadj,solution.poly.nindplus,solution.poly.ngrid,solution.poly.indplus)

    solution.startingguess = false
    solution.alphacoeff = zeros(solution.poly.nfunc*solution.poly.ngrid,2*solution.poly.ns)
    
    #initialize linear solution and kalman matrices
    solution.linsol=initializelinearsolution!(solution.poly.nparams,solution.poly.nvars,solution.poly.nexog,solution.poly.nexogshock,solution.poly.nexogcont,solution.linsol)

    solution.poly.endogsteady=Array{Float64}(undef,solution.poly.nvars+solution.poly.nexog,1)
    solution.poly.slopeconmsv=Array{Float64}(undef,2*nmsvadj,1)
    solution.poly.shockbounds=Array{Float64}(undef,solution.poly.nexogshock,2)
    solution.poly.shockdistance=Array{Float64}(undef,solution.poly.nexogshock,1)
    solution.poly.exoggrid=Array{Float64}(undef,nexogadj,solution.poly.ns)
    
    return 
    
end

function initializetestsolution!(test_solution)
    
    data=Base.DataFmt.readdlm("solution%poly.txt")
    numericData=convert(Array{Int},data[1:end-1,3])
    test_solution.poly.nfunc       =numericData[1]
    test_solution.poly.nmsv        =numericData[2]
    test_solution.poly.nvars       =numericData[3]
    test_solution.poly.ngrid       =numericData[4]
    test_solution.poly.nparams     =numericData[5]
    test_solution.poly.nindplus    =numericData[6]
    test_solution.poly.nexog       =numericData[7]
    test_solution.poly.nexogshock  =numericData[8]
    test_solution.poly.nexogcont   =numericData[9]           
    test_solution.poly.ns          =numericData[10]         
    test_solution.poly.nsexog      =numericData[11]           
    test_solution.poly.ninter      =numericData[12]          
    test_solution.poly.nquad       =numericData[13]
    test_solution.poly.zlbswitch   =true

    data=readdlm("solution%poly%ghnodes.txt")
    test_solution.poly.ghnodes=data

    data=readdlm("solution%poly%ghweights.txt")
    test_solution.poly.ghweights=data'     

    data=readdlm("solution%poly%indplus.txt")
    test_solution.poly.indplus=data

    data=readdlm("solution%poly%nshockgrid.txt",Int)
    test_solution.poly.nshockgrid=data'

    data=readdlm("solution%poly%interpolatemat.txt")
    test_solution.poly.interpolatemat=data #Somthing may be wrong with this ....? Why is it all zeros?

    data=readdlm("solution%poly%endogsteady.txt")
    test_solution.poly.endogsteady=data # I think some of the values are undefined

    data=readdlm("solution%poly%exoggrid.txt")
    test_solution.poly.exoggrid=data 

    data=readdlm("solution%poly%endogsteady.txt")
    test_solution.poly.endogsteady=data 

    data=readdlm("solution%poly%shockbounds.txt")
    test_solution.poly.shockbounds=zeros(size(data)) #something is wrong with the .txt file 

    data=readdlm("solution%poly%shockdistance.txt")
    test_solution.poly.shockdistance=data 

    data=readdlm("solution%poly%slopeconmsv.txt")
    test_solution.poly.slopeconmsv=data 
end
                                       
function gensys(Γ0, Γ1, c, Ψ, Π, args...)
    F = try
        schur!(complex(Γ0), complex(Γ1))
    catch ex
        if isa(ex, LinearAlgebra.LAPACKException)
            info("LAPACK exception thrown while computing Schur decomposition of Γ0 and Γ1.")
            eu = [-3, -3]

            G1 = Array{Float64, 2}(0,0)
            C = Array{Float64, 1}(0)
            impact = Array{Float64, 2}(0,0)
            fmat = Array{Complex{Float64}, 2}(0,0)
            fwt = Array{Complex{Float64}, 2}(0,0)
            ywt = Vector{Complex{Float64}}(0)
            gev = Vector{Complex{Float64}}(0)
            loose = Array{Float64, 2}(0,0)

            return G1, C, impact, fmat, fwt, ywt, gev, eu, loose
        else
            rethrow(ex)
        end
    end
    gensys(F, c, Ψ, Π, args...)
end

function gensys(F::GeneralizedSchur, c, Ψ, Π)
    gensys(F, c, Ψ, Π, new_div(F))
end

# Method that does the real work. Work directly on the decomposition F
function gensys(F::GeneralizedSchur, c, Ψ, Π, div)
    eu = [0, 0]
    ϵ = 1e-6  # small number to check convergence
    nunstab = 0
    zxz = 0
    a, b, = F.S, F.T
    n = size(a, 1)

    select = BitArray(undef,n)
    for i in 1:n
        # nunstab is the variable name used by Chris Sims, but it seems
        # that nunstab should actually correspond to the number of stable λs
        # i.e. nunstab += 1/div > abs(a[i,i])/abs(b[i,i]), which is basically
        # 1 - a small number > abs(a[i,i])/abs(b[i,i])
        select[i] = !(abs(b[i, i]) > div * abs(a[i, i]))
        if (abs(a[i, i]) < ϵ) && (abs(b[i, i]) < ϵ)
            zxz = 1
        end
    end
    nunstab = n - sum(select)

    if zxz == 1
        warn("Coincident zeros. Indeterminacy and/or nonexistence.")
        eu=[-2, -2]

        G1 = Array{Float64, 2}(0, 0)
        C = Array{Float64, 1}(0)
        impact = Array{Float64, 2}(0)
        fmat = Array{Complex{Float64}, 2}(0,0)
        fwt = Array{Complex{Float64}, 2}(0,0)
        ywt = Vector{Complex{Float64}}(0)
        gev = Vector{Complex{Float64}}(0)
        loose = Array{Float64, 2}(0,0)

        return G1, C, impact, fmat, fwt, ywt, gev, eu, loose
    end

    FS = ordschur(F, select)
    a, b, qt, z = FS.S, FS.T, FS.Q, FS.Z
    gev = hcat(diag(a), diag(b))
    qt1 = qt[:, 1:(n - nunstab)]
    qt2 = qt[:, (n - nunstab + 1):n]
    etawt = transpose(conj(qt2))*Π
    neta = size(Π, 2)

    # branch below is to handle case of no stable roots, rather than quitting with an error
    # in that case.
    if nunstab == 0
        etawt = zeros(0, neta)
        ueta = zeros(0, 0)
        deta = zeros(0, 0)
        veta = zeros(neta, 0)
        bigev = 0
    else
        etawtsvd = svd(etawt)
        bigev = findall(etawtsvd.S .> ϵ)
        ueta = etawtsvd.U[:, bigev]
        veta = etawtsvd.V[:, bigev]
        deta = diagm(0 => etawtsvd.S[bigev])
    end

    existence = length(bigev) >= nunstab
    if existence
        eu[1] = 1
    else
        warn("Nonexistence: number of unstable roots exceeds number of jump variables")
    end

    # Note that existence and uniqueness are not just matters of comparing
    # numbers of roots and numbers of endogenous errors.  These counts are
    # reported below because usually they point to the source of the problem.

    # branch below to handle case of no stable roots
    if nunstab == n
        etawt1 = zeros(0, neta)
        bigev = 0
        ueta1 = zeros(0, 0)
        veta1 = zeros(neta, 0)
        deta1 = zeros(0, 0)
    else
        etawt1 = transpose(conj(qt1))*Π
        ndeta1 = min(n - nunstab, neta)
        etawt1svd = svd(etawt1)
        bigev = findall(etawt1svd.S .> ϵ)
        ueta1 = etawt1svd.U[:, bigev]
        veta1 = etawt1svd.V[:, bigev]
        deta1 = diagm(0 => etawt1svd.S[bigev])
    end

    if isempty(veta1)
        unique = true
    else
        loose = veta1 - (veta*transpose(conj(veta))) * veta1
        loosesvd = svd(loose)
        nloose = sum(abs.(loosesvd.S) .> ϵ * n)
        unique = (nloose == 0)
    end

    if unique
        eu[2] = 1
    else
        warn("Indeterminacy: $(nloose) loose endogeneous error(s)")
    end

    tmat = hcat(Matrix{Float64}(I, n - nunstab, n - nunstab), -(ueta * (deta \ veta') * veta1 * 	(deta1*transpose(conj(ueta1))))')

    G0 = vcat(tmat * a, hcat(zeros(nunstab, n - nunstab), Matrix{Float64}(I, nunstab, nunstab)))
    G1 = vcat(tmat * b, zeros(nunstab, n))

    # G0 is always non-singular because by construction there are no zeros on
    # the diagonal of a(1:n-nunstab,1:n-nunstab), which forms G0's ul corner.
    G0I = inv(G0)
    G1 = G0I * G1
    usix = (n - nunstab + 1):n
    Busix = b[usix,usix]
    Ausix = a[usix,usix]
    C = G0I * vcat(tmat * (transpose(conj(qt))*c), (Ausix - Busix) \ (transpose(conj(qt2))* c))
    impact = G0I * vcat(tmat * (transpose(conj(qt))*Ψ), zeros(nunstab, size(Ψ, 2)))
    fmat = Busix \ Ausix
    fwt = -Busix \ (transpose(conj(qt2))* Ψ)
    ywt = G0I[:, usix]

    loose = G0I * vcat(etawt1 * (Matrix{Float64}(I, neta, neta) - (transpose(conj(veta))*veta)), zeros(nunstab, neta))

    G1 = real(z * (G1*transpose(conj(z))))
    C = real(z * C)
    impact = real(z * impact)
    loose = real(z * loose)

    ywt = z * ywt

    return G1, C, impact, fmat, fwt, ywt, gev, eu, loose
end


function new_div(F::GeneralizedSchur)
    ϵ = 1e-6  # small number to check convergence
    n = size(F.T, 1)
    a, b = F.S, F.T
    div = 1.01
    for i in 1:n
        if abs(a[i, i]) > 0
            divhat = abs(b[i, i]) / abs(a[i, i])
            if 1 + ϵ < divhat && divhat <= div
                div = .5 * (1 + divhat)
            end
        end
    end

    return div
end


end
