include("polydef.jl")

using polydef
    #Input user based parameters of the solution
    solution = solutiondetails()
    solution.poly=polydetails() #Must initilize subtype
    solution.linsol=linsoldetails() #Must initilize subtype

    solution.poly.nparams = 43
    solution.poly.nexog = 6
    solution.poly.nexogcont = 0 
    #solution.poly.nexogcont = 1  #adds level technology shock into polynomial 
    solution.poly.nvars = 22
    solution.poly.nmsv = 7
    solution.poly.nfunc = 7
    solution.poly.nindplus = 1

    solution.poly.nshockgrid=[7 3 3 3 3 1]

    npara = solution.poly.nparams
    nvars = solution.poly.nvars+solution.poly.nexog
    nexog = solution.poly.nexog

    if (solution.poly.nindplus == 1) 
        solution.poly.indplus=zeros(Int,solution.poly.nindplus,1)
        solution.poly.indplus[1] = 3
    end

    solution.poly.zlbswitch = true

    #Initilize values for the solution
    initilizesolution!(solution)

    #create a test solution
    test_solution = solutiondetails()
    test_solution.poly=polydetails() #Must initilize subtype
    test_solution.linsol=linsoldetails() #Must initilize subtype

    initilizetestsolution!(test_solution)

    using Base.Test
    @testset "Polydef Tests" begin
        for n in fieldnames(solution.poly)
            tolerance=1e-8
            println("Parameter: ",n)
            @test getfield(solution.poly,n) ≈ getfield(test_solution.poly,n) atol=tolerance
        end
    end
