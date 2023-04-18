using Plots
using LaTeXStrings

#
# Rotina "geral" para fazer os gráficos
#
function Make_plot(tspan::Tuple{Float64,Float64}, dt::Float64, Example::Function, f!::Function, plot_name="plot", beta_c = 1E-2)

    # Generate the discrete times
    tt = tspan[1]:dt:tspan[2]

    # Load Problem data
    M,C,K,U0,V0,t0 = Problem_Data(beta_c)

    # Solve the problem
    y, _ =  Example(M,C,K,U0,V0,t0)

    # Generate the plot
    ndofs = size(y(0.0),1)
    yy = reshape([real(y(t))[k] for k=1:ndofs for t in tt],length(tt),ndofs)

    plot()
    for i=1:ndofs
       plot!(tt, yy[:,i], linewidth = 2, title = "", legend=:outertopright,
           xaxis = L"t", yaxis = L"\mathbf{y}(t)", label=L"y_%$(i)(t)",dpi=500, alpha=0.5)
    end

    # Solve by using Newmark - Beta method
    U,V,A,T = Solve_newmark(M, C, K, f!, tspan, dt, U0=U0, V0=V0)

    #
    # Plot Newmark solution
    #
    for i=1:ndofs
       plot!(T, U[:,i], linewidth = 1, title = "", legend=:outertopright, ls=:dot, 
           xaxis = L"t", yaxis = L"\mathbf{y}(t)", label=L"\tilde{y}_%$(i)(t)",dpi=500, alpha=1.0)
    end


    savefig(plot_name*".pdf")

end

