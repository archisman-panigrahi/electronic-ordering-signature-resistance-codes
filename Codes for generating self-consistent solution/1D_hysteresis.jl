# Define the kinetic energy function and interaction function
# define lattice laplacian function

# Define initial population function
# Write self-consistent method

println("data is already generated in ../data/dataEInducedSwitching.csv")
println("To just see the plot, run plot_E_induced_switching() to display it, since the data is already generated")

println("run generate_data() to regenerate the same data")
println("run plot_E_induced_switching() to generate the plot")
println("Note that it takes ~5 minutes to generate all the data points in a regular desktop computer")


using Plots; using CSV; using DataFrames; using LaTeXStrings
using Plots.PlotMeasures

Plots.default(show = true)

using DelimitedFiles

function kinetic_energy(kx::Float64)
    return 1-2*kx^2+kx^4
end

function kinetic_energy_larger_dist(kx::Float64,dist_kin::Float64)
    return (kx^2-dist_kin^2)^2
end

function Interaction_amp(k1x::Float64, k2x::Float64, V::Float64, q_0::Float64)
    return V*exp(-(k1x-k2x)^2/(2*q_0^2))
end

function lattice_laplacian(f,distgrid)
    lattice_laplacian_f = zero(f)
    for ii = 2:size(f,1)-1
        lattice_laplacian_f[ii] = (f[ii-1] + f[ii+1] - 2*f[ii])/distgrid^2
    end
    return lattice_laplacian_f
end

function OneToLMod(a,L)
    c = mod(a,L)
    if c==0
        return L
    else
        return c
    end
end

function EFieldSwitching(Lx::Int64, V::Float64, mu::Float64, polarized::Float64, Temp::Float64, q_0::Float64, pocket_dist::Float64, E_field::Float64, initial_dist::Vector{Float64}=[1.0])
    ### All the plots are commented out to speed up data generation
    #close("all")
    lattice_points = range(-3.0,3.0,Lx)
    distgrid = lattice_points[2]-lattice_points[1]
    E_k = zeros(Lx)
    lambda_k = zeros(Lx)
    Fermi_occupation = zeros(Lx)
    Fermi_occupation_new = rand(Lx)
    Fermi_occupation_intermediate = rand(Lx)
    Fermi_occupation_shifted = rand(Lx)
    beta = 1/Temp

    shift_int = round(Int64,E_field/distgrid)

    println("V = ", V)
    println("mu = ", mu)

    for ii = 1:Lx
        E_k[ii] = kinetic_energy_larger_dist(lattice_points[ii],pocket_dist)
    end

    ## Instead of this. Try the following: Start with the equilibrium distribution after self-consistency


    if(initial_dist == [1.0])
        if (polarized == 0.0 || polarized == 2.0)
            for ii = 1:Lx 
                Fermi_occupation[ii] = 1/(1+exp(beta*(E_k[ii]-mu)))
            end
        elseif polarized == 1.0
            for ii = Int(round(Lx/2)):Lx
                Fermi_occupation[ii] = 1/(1+exp(beta*(E_k[ii]-mu)))
            end
        elseif polarized == -1.0
            for ii = 1:Int(round(Lx/2))
                Fermi_occupation[ii] = 1/(1+exp(beta*(E_k[ii]-mu)))
            end
        else
            println("Polarization can be 1 (right pocket), -1 (left-pocket), 2 or 0 (two pockets)\n")
            println("You entered polarization = ", polarized, "\n")
            println("Please start over\n")
            close("all")
            return
        end
    else
        Fermi_occupation = copy(initial_dist);
    end

    n_occupied = sum(Fermi_occupation)

    println("n_occupied = ", n_occupied)
    #plt1 = scatter(lattice_points, E_k, title=string("Kinetic energy bandstructure"))
    #plt1 = plot!(lattice_points, mu*ones(Lx))

    #plt2 = scatter(lattice_points, Fermi_occupation, reuse=false, title=string("Initial Occupation function"))

    Fermi_occupation_intermediate = copy(Fermi_occupation)
    Fermi_occupation_shifted = copy(Fermi_occupation)

    this_variable_breaks_while_loop = 0
    iter = 1

    ##Initialize muTrial

    muTrial = 0.0
    while(this_variable_breaks_while_loop < 2)
        lambda_k = zeros(Lx)
        for ii = 1:Lx
            for kk = 1:Lx
                lambda_k[ii] += Interaction_amp(lattice_points[ii], lattice_points[kk], V, q_0) * Fermi_occupation_shifted[kk]
            end
        end

        ##Here we estimate the new Chemical potential
        if iter == 1
            muTrial = copy(mu)
        end

        itermu = 0
        increment = 0.04
        while(true)
            nNew = 0.0
            itermu += 1
            for ii = 1:Lx
                nNew += 1/(1+ exp(beta*(E_k[ii] - lambda_k[ii] - muTrial))) 
            end
            println("nNew = ", nNew, ", itermu = ", itermu, ", n_occupied = ", n_occupied, ", iter = ", iter)

            if(nNew < n_occupied - 0.05)
                muTrial += increment
            elseif(nNew > n_occupied + 0.05)
                muTrial += -increment
            else
                break
            end
            println("iter = ",iter,", itermu = ", itermu, ", muTrial = ", muTrial)

            if(itermu > 25)
                increment = increment/2
                itermu = 0
            end

            if(increment < 10^(-6))
                println("Could not find correct chemical potential")
                break
            end
        end

        #plt3 = scatter(lattice_points, E_k - lambda_k, reuse = false, ylims=(minimum(E_k-lambda_k),1.1), title=string("Intermediate Mean Field bandstructure at iter =", string(iter)))
        #plt3 = plot!(lattice_points, muTrial*ones(Lx))

        for ii=1:Lx
            Fermi_occupation_new[ii] = 1/(1+exp(beta*(E_k[ii] - lambda_k[ii] - muTrial)))
        end

        for ii=1:Lx
            Fermi_occupation_shifted[ii] = Fermi_occupation_new[OneToLMod(ii+shift_int,Lx)]
        end

        #plt4 = scatter(lattice_points, Fermi_occupation_new, reuse=false, title=string("Intermediate Equilibrium Fermi occupation iter = ", string(iter)))
        #plt4pt5 = scatter(lattice_points, Fermi_occupation_shifted, reuse=false, title=string("Intermediate non-Equilibrium Fermi occupation iter = ", string(iter)))

        # Check convergence
        println("shift_int = ", shift_int)
        println("nRight = ", sum(Fermi_occupation_new[Int(round(Lx/2)):Lx]))
        println("nLeft = ", sum(Fermi_occupation_new[1:Int(round(Lx/2))]))
        println("Max error in electron number = ", maximum(abs.(Fermi_occupation_new - Fermi_occupation_intermediate)))
        if maximum(abs.(Fermi_occupation_new - Fermi_occupation_intermediate)) < 0.0008
            this_variable_breaks_while_loop += 1;
        end
        Fermi_occupation_intermediate = copy(Fermi_occupation_new)
        iter += 1

        if(iter > 400)
            println("breaking due to getting stuck")
            break
        end
    end

    #plt5 = scatter(lattice_points, E_k - lambda_k, reuse = false, ylims=(minimum([minimum(E_k-lambda_k),muTrial-0.1]),1.1), title=string("Final Mean Field bandstructure"))
    #plt5 = plot!(lattice_points, muTrial*ones(Lx))

    #plt6 = scatter(lattice_points, Fermi_occupation_new, reuse=false, title=string("Final Mean Field Equ Fermi occupation"))
    #plt6 = scatter(lattice_points, Fermi_occupation_shifted, reuse=false, title=string("Final Mean Field non-Equ Fermi occupation"))


    println("Initial_total_curvature = ",sum(lattice_laplacian(E_k,distgrid).*Fermi_occupation))
    println("Final_total_curvature = ",sum(lattice_laplacian(E_k-lambda_k,distgrid).*Fermi_occupation_new))
    println("------------")
    println("Final_average_curvature = ",sum(lattice_laplacian(E_k-lambda_k,distgrid).*Fermi_occupation_new)/n_occupied)
    println("------------")

    ##Calculation of Free energy
    U = sum((E_k-0.5*lambda_k).*Fermi_occupation_new)
    println("U=",U)
    S = 0

    for ii=1:Lx 
        term = -Fermi_occupation_new[ii] * log(Fermi_occupation_new[ii])
        if term >0
            S += term
        end
        term = -(1-Fermi_occupation_new[ii]) * log(1-Fermi_occupation_new[ii])
        if term >0
            S += term
        end
    end
    Total_E_kinetic = sum(E_k .* Fermi_occupation_new)

    println("TS = ", S/beta)
    println("Total_free_energy = ",U-(S/beta))

    println("------------")
    println("Free_energy_per_particle = ",(U-(S/beta))/(n_occupied))
    println("------------")

    println("'Fermi Energy' from bottom = ", muTrial - minimum(E_k-lambda_k))
    println("Temperature = ",1/beta)

    println("nRight = ", sum(Fermi_occupation_new[Int(round(Lx/2)):Lx]))
    println("nLeft = ", sum(Fermi_occupation_new[1:Int(round(Lx/2))]))
    println("------------")
    println("polarization = ", (sum(Fermi_occupation_new[Int(round(Lx/2)):Lx]) - sum(Fermi_occupation_new[1:Int(round(Lx/2))]))/(sum(Fermi_occupation_new[1:Int(round(Lx/2))]) + sum(Fermi_occupation_new[Int(round(Lx/2)):Lx])))
    println("------------")

    println("iter=",iter)

    #plt5 = scatter(lattice_points,lattice_laplacian(E_k-lambda_k,distgrid), reuse=false, title=string("second derivative"))
    #plt5 = plot!(lattice_points,zeros(Lx))

    #plt6 = scatter(lattice_points,lattice_laplacian(E_k-lambda_k,distgrid).*Fermi_occupation_new, reuse=false, title=string("conductivity"))
    ### Save to file
    dataToSave = [Lx,V,mu,muTrial,polarized,beta,q_0,U,S,U-(S/beta),sum(Fermi_occupation_new[Int(round(Lx/2)):Lx]),sum(Fermi_occupation_new[1:Int(round(Lx/2))]),(sum(Fermi_occupation_new[Int(round(Lx/2)):Lx]) - sum(Fermi_occupation_new[1:Int(round(Lx/2))]))/sum(Fermi_occupation_new),sum(lattice_laplacian(E_k-lambda_k,distgrid).*Fermi_occupation_new),Total_E_kinetic,n_occupied,pocket_dist,E_field]
    io = open("../data/dataEInducedSwitching.csv","a")
    writedlm(io,dataToSave',",")
    close(io)
end

function returnEquilibrium_dist(Lx::Int64, V::Float64, mu::Float64, polarized::Float64, Temp::Float64, q_0::Float64, pocket_dist::Float64)
    lattice_points = range(-3.0,3.0,Lx)
    distgrid = lattice_points[2]-lattice_points[1]
    E_k = zeros(Lx)
    lambda_k = zeros(Lx)
    Fermi_occupation = zeros(Lx)
    Fermi_occupation_new = rand(Lx)
    Fermi_occupation_intermediate = rand(Lx)
    beta = 1/Temp

    println("V = ", V)
    println("mu = ", mu)

    for ii = 1:Lx
        E_k[ii] = kinetic_energy_larger_dist(lattice_points[ii], pocket_dist)
    end

    if (polarized == 0.0 || polarized == 2.0)
        for ii = 1:Lx 
            Fermi_occupation[ii] = 1/(1+exp(beta*(E_k[ii]-mu)))
        end
    elseif polarized == 1.0
        for ii = Int(round(Lx/2)):Lx
            Fermi_occupation[ii] = 1/(1+exp(beta*(E_k[ii]-mu)))
        end
    elseif polarized == -1.0
        for ii = 1:Int(round(Lx/2))
            Fermi_occupation[ii] = 1/(1+exp(beta*(E_k[ii]-mu)))
        end
    else
        println("Polarization can be 1 (one pocket), 2 or 0 (two pockets)\n")
        println("You entered polarization = ", polarized, "\n")
        println("Please start over\n")
        close("all")
        return
    end

    n_occupied = sum(Fermi_occupation)

    println("n_occupied = ", n_occupied)

    Fermi_occupation_intermediate = copy(Fermi_occupation)
    this_variable_breaks_while_loop = 0
    iter = 1

    ##Initialize muTrial

    muTrial = 0.0
    while(this_variable_breaks_while_loop < 2)
        lambda_k = zeros(Lx)
        for ii = 1:Lx
            for kk = 1:Lx
                lambda_k[ii] += Interaction_amp(lattice_points[ii], lattice_points[kk], V, q_0) * Fermi_occupation_intermediate[kk]
            end
        end

        ##Here we estimate the new Chemical potential
        if iter == 1
            muTrial = copy(mu)
        end

        itermu = 0
        increment = 0.04
        while(true)
            nNew = 0.0
            itermu += 1
            for ii = 1:Lx
                nNew += 1/(1+ exp(beta*(E_k[ii] - lambda_k[ii] - muTrial))) 
            end
            println("nNew = ", nNew, ", itermu = ", itermu, ", n_occupied = ", n_occupied, ", iter = ", iter)

            if(nNew < n_occupied - 0.05)
                muTrial += increment
            elseif(nNew > n_occupied + 0.05)
                muTrial += -increment
            else
                break
            end
            println("iter = ",iter,", itermu = ", itermu, ", muTrial = ", muTrial)

            if(itermu > 25)
                increment = increment/2
                itermu = 0
            end

            if(increment < 10^(-6))
                println("Could not find correct chemical potential")
                break
            end
        end

        for ii=1:Lx
            Fermi_occupation_new[ii] = 1/(1+exp(beta*(E_k[ii] - lambda_k[ii] - muTrial)))
        end

        # Check convergence
        println("nRight = ", sum(Fermi_occupation_new[Int(round(Lx/2)):Lx]))
        println("nLeft = ", sum(Fermi_occupation_new[1:Int(round(Lx/2))]))
        println("Max error in electron number = ", maximum(abs.(Fermi_occupation_new - Fermi_occupation_intermediate)))
        if maximum(abs.(Fermi_occupation_new - Fermi_occupation_intermediate)) < 0.0008
            this_variable_breaks_while_loop += 1;
        end
        Fermi_occupation_intermediate = copy(Fermi_occupation_new)
        iter += 1

        if(iter > 400)
            println("breaking due to getting stuck")
            break
        end
    end
    return Fermi_occupation_new
end


### First run with small n to test for errors
function generate_data()
    n=101
    E_field_array = range(-0.7,0.7,n)

    f_initial_array_left = returnEquilibrium_dist(2001,0.01,0.352,-1.0,0.2,0.4,1.0)
    f_initial_array_right = returnEquilibrium_dist(2001,0.01,0.352,1.0,0.2,0.4,1.0)

    for ii = 1:n
        EFieldSwitching(2001,0.01,0.352,-1.0,0.2,0.4,1.0,E_field_array[ii],f_initial_array_left)
    end

    for ii = 1:n
        EFieldSwitching(2001,0.01,0.352,1.0,0.2,0.4,1.0,E_field_array[ii],f_initial_array_right)
    end
end

function plot_E_induced_switching()
    n=101
    data = CSV.read("../data/dataEInducedSwitching.csv",DataFrame,header=0);
    polarization_array_Left = data[1:n,end-5];
    polarization_array_Right = data[n+1:2*n,end-5];

    E_Field_array_Left = data[1:n,end];
    E_Field_array_Right = data[n+1:2*n,end];

    plt1 = scatter(E_Field_array_Right,polarization_array_Right,framestyle = :box,legend=false, grid=false)
    plt1 = scatter!(E_Field_array_Left,polarization_array_Left,framestyle = :box,legend=false, grid=false)


    gr();
    plt1 = scatter(E_Field_array_Right,polarization_array_Right,framestyle = :box,legend=false, grid=false,
    marker = :circle,
    markercolor = :red,
    markerstrokecolor = :red,
    markersize=4,
    size = (800, 600),
    xlabel=L"e \tau E/\hbar k_0",
    ylabel=string(L"\textrm{Pocket}"," ",L"\textrm{polarization}"," ",L"(\phi_{1D})"),
    xguidefontsize=20,
    yguidefontsize=18,
    xtickfont = font(18),
    ytickfont = font(18),
    bottom_margin=6.5px
    )
    plt1=xticks!(-0.5:0.25:0.5, [L"$-0.50$", L"$-0.25$", L"$0.00$", L"$0.25$", L"$0.50$"])
    plt1=yticks!(-1:0.5:1.0, [L"$-1.0$",L"$-0.5$",L"$0.0$", L"$0.5$", L"$1.0$"])
    plt1 = plot!(E_Field_array_Right,polarization_array_Right, linestyle=:dash, linecolor=:red, linewidth=:1.0)


    plt1 = scatter!(E_Field_array_Left,polarization_array_Left,framestyle = :box,legend=false, grid=false,
    marker = :circle,
    markercolor = :blue,
    markerstrokecolor = :blue,
    markersize=2
    )
    plt1=xticks!(-0.5:0.25:0.5, [L"$-0.50$", L"$-0.25$", L"$0.00$", L"$0.25$", L"$0.50$"])
    plt1=yticks!(-1:0.5:1.0, [L"$-1.0$",L"$-0.5$",L"$0.0$", L"$0.5$", L"$1.0$"])
    plt1 = plot!(E_Field_array_Left,polarization_array_Left, linestyle=:dash, linecolor=:blue, linewidth=:1.0)

    savefig("E_induced_switching.pdf")
    display(plt1)
end