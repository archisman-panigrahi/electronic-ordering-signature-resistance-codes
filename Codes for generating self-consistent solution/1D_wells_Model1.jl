println("How to run:\n mainfunction(Lx, V, mu, polarized, Temp, q_0)")
println("For example")
println("@time mainfunction(1001,0.01,0.2,1.0,0.1,0.4)")

using Plots
Plots.default(show = true)
pyplot()

using DelimitedFiles

function kinetic_energy(kx::Float64)
    return 1-2*kx^2+kx^4
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


function mainfunction(Lx::Int64, V::Float64, mu::Float64, polarized::Float64, Temp::Float64, q_0::Float64)
    close("all")
    lattice_points = range(-1.5,1.5,Lx)
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
        E_k[ii] = kinetic_energy(lattice_points[ii])
    end

    if (polarized == 0.0 || polarized == 2.0)
        for ii = 1:Lx 
            Fermi_occupation[ii] = 1/(1+exp(beta*(E_k[ii]-mu)))
        end
    elseif polarized == 1.0
        for ii = Int(round(Lx/2)):Lx
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
    plt1 = scatter(lattice_points, E_k, title=string("Kinetic energy bandstructure"))
    plt1 = plot!(lattice_points, mu*ones(Lx))

    plt2 = scatter(lattice_points, Fermi_occupation, reuse=false, title=string("Initial Occupation function"))

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

        plt3 = scatter(lattice_points, E_k - lambda_k, reuse = false, title=string("Intermediate Mean Field bandstructure at iter =", string(iter)))
        plt3 = plot!(lattice_points, muTrial*ones(Lx))

        for ii=1:Lx
            Fermi_occupation_new[ii] = 1/(1+exp(beta*(E_k[ii] - lambda_k[ii] - muTrial)))
        end

        plt4 = scatter(lattice_points, Fermi_occupation_new, reuse=false, title=string("Intermediate Mean Field Fermi occupation at iter = ", string(iter)))

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

    plt3 = scatter(lattice_points, E_k - lambda_k, reuse = false, title=string("Final Mean Field bandstructure"))
    plt3 = plot!(lattice_points, muTrial*ones(Lx))

    plt4 = scatter(lattice_points, Fermi_occupation_new, reuse=false, title=string("Final Mean Field Fermi occupation"))

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

    plt5 = scatter(lattice_points,lattice_laplacian(E_k-lambda_k,distgrid), reuse=false, title=string("second derivative"))
    plt5 = plot!(lattice_points,zeros(Lx))

    plt6 = scatter(lattice_points,lattice_laplacian(E_k-lambda_k,distgrid).*Fermi_occupation_new, reuse=false, title=string("conductivity"))
    ### Save to file
    #[Lx,V,muInitial,muFinal,polarized,beta,q_0,U,S,Free_Energy,nRight,nLeft,(nRight-nLeft)/nTotal,Final_conductivity,total_kinetic,nTotal]
    dataToSave = [Lx,V,mu,muTrial,polarized,beta,q_0,U,S,U-(S/beta),sum(Fermi_occupation_new[Int(round(Lx/2)):Lx]),sum(Fermi_occupation_new[1:Int(round(Lx/2))]),(sum(Fermi_occupation_new[Int(round(Lx/2)):Lx]) - sum(Fermi_occupation_new[1:Int(round(Lx/2))]))/sum(Fermi_occupation_new),sum(lattice_laplacian(E_k-lambda_k,distgrid).*Fermi_occupation_new),Total_E_kinetic,n_occupied]
    io = open("dataTvsN2pocketT=0.0142.csv","a")
    writedlm(io,dataToSave',",")
    close(io)
end
