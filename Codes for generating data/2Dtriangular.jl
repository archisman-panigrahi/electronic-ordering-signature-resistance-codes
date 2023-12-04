###Example
### @time mainfunction(10,5,2.0,0.4,-2.0,2.0,50.0,1.5,"data2DTriangular_1pocket_initial.csv")

using Plots
Plots.default(show = true)
pyplot()

using DelimitedFiles

include("triangular_lattice_generator.jl")

function kinetic_energy(kx::Float64,ky::Float64,pocket1_x::Float64,pocket1_y::Float64,pocket2_x::Float64,pocket2_y::Float64,pocket3_x::Float64,pocket3_y::Float64,depth_of_pocket::Float64,radius_of_pocket::Float64)
    return -depth_of_pocket*(exp(-((kx-pocket1_x)^2 + (ky-pocket1_y)^2)/(2*radius_of_pocket^2)) + exp(-((kx-pocket2_x)^2 + (ky-pocket2_y)^2)/(2*radius_of_pocket^2)) + exp(-((kx-pocket3_x)^2 + (ky-pocket3_y)^2)/(2*radius_of_pocket^2)))
end

function Interaction_amp(k1x::Float64, k1y::Float64, k2x::Float64, k2y::Float64, V::Float64, q_0::Float64)
    return V*exp(-((k1x-k2x)^2 + (k1y-k2y)^2)/(2*q_0^2))
end

function lattice_laplacian_triangular(f::Vector{Float64},coordinates_x::Vector{Float64},coordinates_y::Vector{Float64},distgrid::Float64)
    lattice_laplacian_f = zero(f) # initialize to zero since we don't touch the boundaries
    n_neighbors = zero(f)
    sum = zero(f)
    distance = 0.0
    for ii = 1:size(f,1)
        for jj = 1:size(f,1)
            distance = sqrt((coordinates_x[ii]-coordinates_x[jj])^2 + (coordinates_y[ii]-coordinates_y[jj])^2)
            #println("The distance between ",ii,"and ",jj," is", distance)
            if(distance > 0 && distance < 1.1*distgrid) ##We do this to get rid of floating point round errors
                sum[ii]+=f[jj]
                n_neighbors[ii] += 1
                #println("The neighbor of ",ii,"is ",jj)
            end
        end
        if(n_neighbors[ii] == 6)
            lattice_laplacian_f[ii] = (2/3)*(sum[ii] - 6*f[ii])/distgrid^2
        end
    end
    return lattice_laplacian_f
end

function check_lattice_laplacian(n::Int64,a::Float64)
    close("all")
    nSites = 3*n^2 - 3*n + 1
    coordinates_x, coordinates_y = generate_triangular_lattice_on_a_hexagon(n,a);

    plt1=scatter(coordinates_x,coordinates_y, legend=false, aspect_ratio=:equal)

    fFunction = coordinates_x.^2 + coordinates_y.^2;

    plt2 = scatter(coordinates_x,coordinates_y, zcolor=fFunction, markersize=10, colormap=:rainbow, reuse=false,aspect_ratio=:equal)
    lat_lap_f = lattice_laplacian_triangular(fFunction,coordinates_x,coordinates_y,a);
    plt3 = scatter(coordinates_x,coordinates_y, zcolor=lat_lap_f, markersize=10, colormap=:rainbow, reuse=false,aspect_ratio=:equal, title="Numerical Laplacian")
    plt4 = scatter(coordinates_x,coordinates_y, lat_lap_f, markersize=10, reuse=false,aspect_ratio=:equal, title="Numerical Laplacian")
end

function second_derivative_x(f::Vector{Float64},coordinates_x::Vector{Float64},coordinates_y::Vector{Float64},distgrid::Float64)
    x_second_derivative = zero(f) # initialize to zero since we don't touch the boundaries
    n_neighbors = zero(f)
    sum = zero(f)
    distance_x = 0.0
    distance_y = 0.0

    min_y = minimum(coordinates_y)
    max_y = maximum(coordinates_y)

    for ii = 1:size(f,1)
        for jj = 1:size(f,1)
            distance_x = abs(coordinates_x[ii]-coordinates_x[jj])
            distance_y = abs(coordinates_y[ii]-coordinates_y[jj])
            #println("The distance between ",ii,"and ",jj," is", distance)
            if(coordinates_y[ii] > min_y && coordinates_y[ii] < max_y && distance_x > 0 && distance_x < 1.1*distgrid && distance_y < 0.1*distgrid) ##We do this to get rid of floating point round errors
                sum[ii]+=f[jj]
                n_neighbors[ii] += 1
                #println("The neighbor of ",ii,"is ",jj)
            end
        end
        if(n_neighbors[ii] == 2)
            x_second_derivative[ii] = (sum[ii] - 2*f[ii])/distgrid^2
        end
    end
    return x_second_derivative
end

function check_x_second_derivative(n::Int64,a::Float64)
    close("all")
    nSites = 3*n^2 - 3*n + 1
    coordinates_x, coordinates_y = generate_triangular_lattice_on_a_hexagon(n,a);

    plt1=scatter(coordinates_x,coordinates_y, legend=false, aspect_ratio=:equal)

    fFunction = exp.(-((coordinates_x - sum(coordinates_x) * ones(nSites)/nSites).^2 +(coordinates_y - sum(coordinates_y)* ones(nSites)/nSites).^2)/3)

    plt2 = scatter(coordinates_x,coordinates_y, fFunction, markersize=10, reuse=false,aspect_ratio=:equal)
    x_second_derivative = second_derivative_x(fFunction,coordinates_x,coordinates_y,a);
    plt3 = scatter(coordinates_x,coordinates_y, x_second_derivative, markersize=10, reuse=false,aspect_ratio=:equal, title="Numerical second derivative")
    plt3 = scatter!(coordinates_x,coordinates_y, ((4/9)*(coordinates_x - sum(coordinates_x) * ones(nSites)/nSites).^2 - (2/3)*ones(nSites)) .* exp.(-((coordinates_x - sum(coordinates_x) * ones(nSites)/nSites).^2 +(coordinates_y - sum(coordinates_y)* ones(nSites)/nSites).^2)/3), markersize=6, reuse=true,aspect_ratio=:equal)
end

function dipole_moment(n::Int64, fermi_function::Vector{Float64},coordinates_x::Vector{Float64},coordinates_y::Vector{Float64})
    nSites =  3*n^2 - 3*n + 1
    center_index = Int((nSites - (2*n-1))/2 + n)
    distgrid = coordinates_x[2] - coordinates_x[1]
    dipole_moment_x = 0.0
    dipole_moment_y = 0.0

    for ii = 1:nSites
        dipole_moment_x += distgrid * (coordinates_x[ii]-coordinates_x[center_index]) * fermi_function[ii]
        dipole_moment_y += distgrid * (coordinates_y[ii]-coordinates_y[center_index]) * fermi_function[ii]
    end

    return dipole_moment_x, dipole_moment_y
end

function mainfunction(n::Int64, dist_pockets_from_center::Int64, rad_pockets::Float64,V::Float64, mu::Float64, polarized::Float64, beta::Float64, q_0::Float64, outputfilename::AbstractString)
    close("all")
    nSites = 3*n^2 - 3*n + 1
    coordinates_x, coordinates_y = generate_triangular_lattice_on_a_hexagon(n,1.0)
    distgrid = 1.0
    E_k = zeros(nSites)
    lambda_k = zeros(nSites)
    Fermi_occupation = zeros(nSites)
    Fermi_occupation_new = rand(nSites)
    Fermi_occupation_intermediate = rand(nSites)

    center_index = Int((nSites - (2*n-1))/2 + n)
    triangular_1_index = Int(center_index -dist_pockets_from_center)
    triangular_2_index = Int(n + n*(n-dist_pockets_from_center-1) + (n-dist_pockets_from_center-2)*(n-dist_pockets_from_center-1)/2)
    triangular_3_index = Int(center_index + (2*n-1) + (2*n-1)*(dist_pockets_from_center-1) - dist_pockets_from_center*(dist_pockets_from_center-1)/2)

    array1 = [coordinates_x[triangular_1_index],coordinates_x[triangular_2_index],coordinates_x[triangular_3_index]]
    array2 = [coordinates_y[triangular_1_index],coordinates_y[triangular_2_index],coordinates_y[triangular_3_index]]

    plt1=scatter(coordinates_x,coordinates_y, legend=false, aspect_ratio=:equal)
    plt1=scatter!(array1,array2, markersize=8, color=:red)
    plt1=scatter!([coordinates_x[center_index]],[coordinates_y[center_index]], markersize=7, color=:green)
    display(#plt1)
    conductivity_at_first_iteration = 0;
    for ii = 1:nSites
        E_k[ii] = kinetic_energy(coordinates_x[ii],coordinates_y[ii],coordinates_x[triangular_1_index],coordinates_y[triangular_1_index],coordinates_x[triangular_2_index],coordinates_y[triangular_2_index],coordinates_x[triangular_3_index],coordinates_y[triangular_3_index],1.0,rad_pockets)
    end

    if (polarized == 0.0 || polarized == 3.0)
        for ii = 1:nSites
            Fermi_occupation[ii] = 1/(1+exp(beta*(E_k[ii]-mu)))
        end
    elseif polarized == 1.0
        for ii = 1:nSites
            if(coordinates_y[ii] >= sqrt(3)*(coordinates_x[ii]-coordinates_x[1]) -0.01*distgrid && coordinates_y[ii] <= coordinates_y[center_index]-sqrt(3)*(coordinates_x[ii]-coordinates_x[center_index]) + 0.01*distgrid)
                Fermi_occupation[ii] = 1/(1+exp(beta*(E_k[ii]-mu)))
            end
        end
    elseif polarized == 2.0
        for ii = 1:nSites
            Fermi_occupation[ii] = 1/(1+exp(beta*(E_k[ii]-mu)))
        end
        for ii = 1:nSites
            if(coordinates_y[ii] >= sqrt(3)*(coordinates_x[ii]-coordinates_x[1]) -0.01*distgrid && coordinates_y[ii] <= coordinates_y[center_index]-sqrt(3)*(coordinates_x[ii]-coordinates_x[center_index]) + 0.01*distgrid)
                Fermi_occupation[ii] = 0.0
            end
        end
    else
        println("Polarization can be 1 (one pocket), 2 (two pockets), 3 or 0 (three pockets, i.e. unpolarized)\n")
        println("You entered polarization = ", polarized, "\n")
        println("Please start over\n")
        close("all")
        return
    end

    n_occupied = sum(Fermi_occupation)
    println("n_occupied = ", n_occupied)

    plt2 = scatter(coordinates_x,coordinates_y, zcolor=E_k, markersize=8, reuse=false, aspect_ratio=:equal, title=string("Kinetic energy bandstructure"))

    plt2pt5 = scatter(coordinates_x,coordinates_y, E_k, markersize=8, reuse=false, aspect_ratio=:equal, title=string("Kinetic energy bandstructure"))
    plt2pt5 = scatter!(coordinates_x,coordinates_y, mu*ones(nSites), markersize=1, reuse=true, aspect_ratio=:equal, title=string("Kinetic energy bandstructure"))


    plt3 = scatter(coordinates_x,coordinates_y, zcolor=Fermi_occupation, markersize=8, reuse=false, aspect_ratio=:equal, title=string("Initial Occupation function"))

    plt3pt5 = scatter(coordinates_x,coordinates_y, Fermi_occupation, markersize=8, reuse=false, aspect_ratio=:equal, title=string("Initial Occupation function"))


    Fermi_occupation_intermediate = copy(Fermi_occupation)
    this_variable_breaks_while_loop = 0
    iter = 1

    muTrial = copy(mu)

    while(this_variable_breaks_while_loop<2)
        lambda_k = zeros(nSites)
        for ii = 1:nSites
            for jj = 1:nSites
                lambda_k[ii] += Interaction_amp(coordinates_x[ii],coordinates_y[ii],coordinates_x[jj],coordinates_y[jj],V,q_0) * Fermi_occupation_intermediate[jj]
            end
        end

        ##Find the new chemical potential

        itermu = 0
        increment = 0.5

        while(true)
            nNew = 0.0
            itermu += 1
            for ii=1:nSites
                nNew += 1/(1+exp(beta*(E_k[ii]-lambda_k[ii]-muTrial)))
            end
            #If convergence is not tight here, then the Fermi function does not converge easily as the electron number keeps fluctuating
            if(nNew < n_occupied - 0.025)
                muTrial += increment
            elseif(nNew > n_occupied + 0.025)
                muTrial += -increment
            else
                break
            end

            println("iter = ", iter, ", itermu = ", itermu, ", muTrial = ", muTrial)
            println("error in electron number = ", abs(nNew-n_occupied))

            if(itermu > 20)
                increment = increment/2
                itermu = 0
            end

            if(increment < 10^(-6))
                println("Could not find correct chemical potential")
                break
            end
        end

        for ii = 1:nSites
            Fermi_occupation_new[ii] = 1/(1+exp(beta*(E_k[ii]-lambda_k[ii]-muTrial)))
        end

        if(iter == 1)
            #pltiter1 = scatter(coordinates_x,coordinates_y, (E_k - lambda_k), markersize=6, reuse=false, aspect_ratio=:equal, title=string("iter=1 bandstructure"))
            #pltiter1 = scatter!(coordinates_x,coordinates_y, muTrial*ones(nSites), markersize=1, reuse=true, aspect_ratio=:equal, title=string("iter=1 bandstructure"))
            conductivity_at_first_iteration = sum(lattice_laplacian_triangular(E_k-lambda_k,coordinates_x,coordinates_y,1.0).*Fermi_occupation_new)
            #pltdummy = scatter(coordinates_x,coordinates_y, reuse=false, aspect_ratio=:equal, title=string("dummy plot to be replaced"))

        end
        #plt4 = scatter(coordinates_x,coordinates_y, zcolor=Fermi_occupation_new, markersize=8, reuse=true, aspect_ratio=:equal, title=string("Intermediate Mean Field Fermi occupation at iter = ", string(iter)))

        ###Code for printing electron number in the pockets goes here
        #####
        #####

        if maximum(abs.(Fermi_occupation_new - Fermi_occupation_intermediate)) < 0.0005
            this_variable_breaks_while_loop += 1;
        end
        println("iter = ", iter, "max change in Fermi function = ", maximum(abs.(Fermi_occupation_new - Fermi_occupation_intermediate)))
        Fermi_occupation_intermediate = copy(Fermi_occupation_new)
        iter += 1

        if(iter > 400)
            println("Breaking due to getting stuck")
            break
        end
    end

    plt6 = scatter(coordinates_x,coordinates_y, zcolor=(E_k - lambda_k), markersize=8, reuse=false, aspect_ratio=:equal, title=string("final bandstructure"))
    plt6pt5 = scatter(coordinates_x,coordinates_y, (E_k - lambda_k), markersize=6, reuse=false, aspect_ratio=:equal, title=string("final bandstructure"))
    plt6pt5 = scatter!(coordinates_x,coordinates_y, muTrial*ones(nSites), markersize=0.5, reuse=true, aspect_ratio=:equal, title=string("final bandstructure"))


    plt7 = scatter(coordinates_x,coordinates_y, zcolor=Fermi_occupation_new, markersize=8, reuse=false, aspect_ratio=:equal, title=string("Final Occupation function"))

    println("Initial_conductivity = ", sum(lattice_laplacian_triangular(E_k,coordinates_x,coordinates_y,1.0).*Fermi_occupation))
    println("1st iteration conductivity = ", conductivity_at_first_iteration)
    final_conductivity = sum(lattice_laplacian_triangular(E_k-lambda_k,coordinates_x,coordinates_y,1.0).*Fermi_occupation_new)

    ##longitudinal_conductivity sigma_xx is parallel to the polarization direction
    longitudinal_conductivity = sum(second_derivative_x(E_k-lambda_k,coordinates_x,coordinates_y,1.0).*Fermi_occupation_new)
    ##transeverse conductivity sigma_yy is perpendicular to the polarization direction
    transverse_conductivity = final_conductivity - longitudinal_conductivity

    println("Final_laplacian_conductivity = ", final_conductivity)
    println("-----------")
    println("Final_average_curvature = ", final_conductivity/sum(Fermi_occupation_new))
    println("-----------")

    println("Longitudinal conductivity = ", longitudinal_conductivity)
    println("Transverse conductivity = ", transverse_conductivity)

    U = sum((E_k-0.5*lambda_k).*Fermi_occupation_new)
    println("U=",U)
    S = 0
    for ii=1:nSites
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
    println("-----------")
    println("Free energy per electron = ", (U-(S/beta))/sum(Fermi_occupation_new))
    println("-----------")

    println("'Fermi Energy' from bottom = ", muTrial - minimum(E_k-lambda_k))
    println("Final chemical potential = ", muTrial)
    println("-----------")
    println("Temperature = ",1/beta)
    println("-----------")

    println("Initial number of electrons = ", n_occupied)
    println("-----------")
    println("Final number of electrons = ", sum(Fermi_occupation_new))
    println("-----------")
    println("iter=",iter)
    dipole_moment_x, dipole_moment_y = dipole_moment(n,Fermi_occupation_new, coordinates_x, coordinates_y)
    scaled_dipole_moment = sqrt(dipole_moment_x^2 + dipole_moment_y^2)/(sum(Fermi_occupation_new)*distgrid*sqrt((coordinates_x[triangular_1_index]-coordinates_x[center_index])^2 + (coordinates_y[triangular_1_index]-coordinates_y[center_index])^2))
    println("-----------")
    println("scaled dipole moment = ", scaled_dipole_moment)
    println("-----------")

    println("exchange/kinetic = ", abs((U-Total_E_kinetic)/Total_E_kinetic))
    ### Save to file
    #[n,dist_pockets_from_center,rad_pockets,V,muInit,mu_Final,polarized,beta,q_0,Total_energy,S,Free_energy,Total_Energy_kinetic, total_electrons, final_conductivity, scaled_dipole_moment]
    dataToSave = [n,dist_pockets_from_center,rad_pockets,V,mu,muTrial,polarized,beta,q_0,U,S,U-(S/beta),Total_E_kinetic, sum(Fermi_occupation_new), final_conductivity, scaled_dipole_moment, longitudinal_conductivity, transverse_conductivity]
    io = open(outputfilename,"a")
    writedlm(io,dataToSave',",")
    close(io)
end
