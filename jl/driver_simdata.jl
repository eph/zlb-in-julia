include("src/model/class_model.jl")
    import Main.polydef: polydetails, linsoldetails, solutiondetails
    import Main.class_model: model, new_model, solve_serial, describe!, simulate_modeldata
    import Main.get_decisionrule: nonlinearsolver
    
    using LinearAlgebra
    using DelimitedFiles

    #tic()

    #Initilize variables
    parallelswitch = false
    dum = zeros(1)

    #if (parallelswitch .eqv. .true.) then
       #call MPI_init(mpierror)
       #call MPI_Comm_size(MPI_COMM_WORLD,nproc,mpierror)
       #call MPI_Comm_rank(MPI_COMM_WORLD,rank,mpierror)
       #write(*,*) 'Hello from processor ', rank, 'I am 1 of ', nproc ,' processes running'
    #else
       rank = 0
       println("You are running the serial version of the code.")
    #end if
    
    ndsets=1
    capt_simdata=1000 #Should be set to 1000
    modelindex=1

    #iniitialize solution details
    if (modelindex == 0)
        zlbswitch = true 
    else
        zlbswitch = false
    end
    m = new_model(zlbswitch)  
    if (rank == 0)
        describe!(m)
        println("capt = ", capt_simdata)
        println("number of dsets = ", ndsets)
        println("----------------------------")
    else
	println("rank = ", rank)
    end

    #allocate matrices for zlb frequency and duration data
    modeldata=Array{Float64}(undef,m.solution.poly.nvars+2*m.solution.poly.nexog,capt_simdata)

  
    ndsets_used = 0
    for id in 1:ndsets

        #get parameters from disk
        #write(filename,"(A,I4.4,A)") 'results/thinned_posterior/parasim', id-1 , '.txt'
        filename = "input/mean.txt"
        m.params=readdlm(filename,Float64)
 
        #solve model
        #if (parallelswitch == true)
            #convergence = m.solve_parallel(m.params,nproc,rank)
        #else    
            m.solution, convergence = solve_serial(m,m.params) #ACTUAL CODE THAT SOVLES FOR ALPHA PARAMATERS
        #end

        if (convergence  == false) #if no solution, report this back to disk
            if (rank == 0)
                println("Failed to converge for parameters when id = ", id)
                println("----------------------------------")
                if (modelindex == 0)
                    filename="results/zlbstat/modeldata_file"*string(id-1)*".txt"
                else
                    filename="results/zlbstat/modeldata_unc_file"*string(id-1)*".txt"
                end
                writedlm(filename,dum)
                if (modelindex == 0)
                    filename = "results/zlbstat/modeldata_linear_file"*string(id-1)*".txt"
                    writedlm(filename,dum)
                end
x            end
        else  #if computed solution, simulate and send zlb data to disk
            if (rank == 0)
                println("Successfully solved model when id = ", id) 
                println("---------------------------------")
                global ndsets_used = ndsets_used + 1
                nonlinearswitch = true
                seed=1221+id
                modeldata = simulate_modeldata(m,capt_simdata,nonlinearswitch,seed)
                if (modelindex == 0)
                    filename="results/zlbstat/modeldata_file"*string(id-1)*".txt"
                else
                    filename="results/zlbstat/modeldata_unc_file"*string(id-1)*".txt"
                end
                writedlm(filename,modeldata)  
                if (modelindex == 0)
                    nonlinearswitch = false
                    modeldata = simulate_modeldata(m,capt_simdata,nonlinearswitch,seed)                    
                    filename="results/zlbstat/modeldata_linear_file"*string(id-1)*".txt"
                    writedlm(filename,modeldata) #Write matrix array to filename
                end
            end
        end
    end

    if (rank == 0) 
        println("Number of dsets for selected moments table = ", ndsets_used)  
    end
  

#toc()
