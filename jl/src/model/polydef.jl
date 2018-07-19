module polydef
    import FastGaussQuadrature: gausshermite
    export linsoldetails, polydetails, solutiondetails, initializesolution!, setgridsize, 
        exoggridindex, ghquadrature,sparsegrid, smolyakpoly, initializelinearsolution!, initializetestsolution!

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
                exoggridindex[ie,ngrid[ie]*blocksize*(ic-1)+blocksize*(ib-1)+1:ngrid[ie]*blocksize*(ic-1)+blocksize*ib] = ib
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
    ghweights=Array{Float64}(nquadsingle^nexog,1)#Not sure I should make this zeros
    const const_pi = 3.14159265358979323846
    
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
                ghnodes[ie,left:right] = sqrt(2)*quadnodes_s[ib]
                ghweights_mat[ie,left:right] = quadweights_s[ib]
            end
        end
    end

    ghweights = (1.0/const_pi)^(nexog/2.0)*prod(ghweights_mat,1) # product of ghweights_mat along the first dimension

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
    smolyakpoly=Array{Float64}(ngrid,1)

    smolyakpoly[1] = 1.0
    for i in 1:nmsv
        smolyakpoly[2*i] = xx[i]
        smolyakpoly[2*i+1] = 2.0*xx[i]^2-1.0
    end

    for i in 1:nindplus
        smolyakpoly[2*nmsv+2*(i-1)+2] = 4.0*xx[indplus[i]]^3-3.0*xx[indplus[i]]
        smolyakpoly[2*nmsv+2*(i-1)+3] = 8.0*xx[indplus[i]]^4-8.0*xx[indplus[i]]^2+1.0
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
        bbt[:,i] = smolyakpoly(nmsv,ngrid,nindplus,indplus,xgrid[:,i])
    end

    
    # find bbt inverse matrix
    bbtinv = copy(bbt)
    bbtinv,ipiv,info=Base.LinAlg.LAPACK.getrf!(bbtinv)
    if (info == 0) 
        Base.LinAlg.LAPACK.getri!(bbtinv,ipiv)
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
    
    linsol.pp=Array{Float64}(linsol.nvars,linsol.nvars)
    linsol.sigma=Array{Float64}(linsol.nvars,linsol.nexog)
    linsol.endogsteady=Array{Float64}(linsol.nvars,1)
        
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
    solution.exogvarinfo=Array{Int64}(nexogadj,solution.poly.ns)
    solution.exogvarinfo[1:solution.poly.nexogshock,:] = exoggridindex(solution.poly.nshockgrid,solution.poly.nexogshock,solution.poly.ns)
    solution.exogvarinfo[solution.poly.nexogshock+1:nexogadj,:] = 1

    #get matrix used for interpolating the shocks
    solution.poly.interpolatemat=Array{Int64}(solution.poly.nexogshock,2^solution.poly.nexogshock)
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
                solution.poly.interpolatemat[i,left:right] =k-1 
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

    solution.poly.endogsteady=Array{Float64}(solution.poly.nvars+solution.poly.nexog,1)
    solution.poly.slopeconmsv=Array{Float64}(2*nmsvadj,1)
    solution.poly.shockbounds=Array{Float64}(solution.poly.nexogshock,2)
    solution.poly.shockdistance=Array{Float64}(solution.poly.nexogshock,1)
    solution.poly.exoggrid=Array{Float64}(nexogadj,solution.poly.ns)
    
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

    data=Base.DataFmt.readdlm("solution%poly%ghnodes.txt")
    test_solution.poly.ghnodes=data

    data=Base.DataFmt.readdlm("solution%poly%ghweights.txt")
    test_solution.poly.ghweights=data'     

    data=Base.DataFmt.readdlm("solution%poly%indplus.txt")
    test_solution.poly.indplus=data

    data=Base.DataFmt.readdlm("solution%poly%nshockgrid.txt",Int)
    test_solution.poly.nshockgrid=data'

    data=Base.DataFmt.readdlm("solution%poly%interpolatemat.txt")
    test_solution.poly.interpolatemat=data #Somthing may be wrong with this ....? Why is it all zeros?

    data=Base.DataFmt.readdlm("solution%poly%endogsteady.txt")
    test_solution.poly.endogsteady=data # I think some of the values are undefined

    data=Base.DataFmt.readdlm("solution%poly%exoggrid.txt")
    test_solution.poly.exoggrid=data 

    data=Base.DataFmt.readdlm("solution%poly%endogsteady.txt")
    test_solution.poly.endogsteady=data 

    data=Base.DataFmt.readdlm("solution%poly%shockbounds.txt")
    test_solution.poly.shockbounds=zeros(size(data)) #something is wrong with the .txt file 

    data=Base.DataFmt.readdlm("solution%poly%shockdistance.txt")
    test_solution.poly.shockdistance=data 

    data=Base.DataFmt.readdlm("solution%poly%slopeconmsv.txt")
    test_solution.poly.slopeconmsv=data 
end


end
