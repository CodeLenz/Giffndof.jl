using Plots
using LaTeXStrings
#
# Create the plots like in the manuscript
#
function Make_plot(T::Vector{T1},Y::AbstractMatrix{T1},U::AbstractMatrix{T1},plot_name="plot";dofs=[]) where T1

     
    # Dimensions
    ndofs = size(U,1)

    # Create the plot
    plot()

    # If dofs is empty, plot all dofs
    if isempty(dofs)
        dofs = 1:ndofs
    end

    # Plot Y
    for i in dofs
        plot!(T, Y[:,i], linewidth = 2, title = "", legend=:outertopright,
            xaxis = L"t", yaxis = L"\mathbf{y}(t)", label=L"y_%$(i)(t)",dpi=500, alpha=1.0)
    end

    #
    # Plot Newmark solution
    #
    for i in dofs
        plot!(T, U[i,:], linewidth = 1, title = "", legend=:outertopright, ls=:dot, palette=:thermal, 
            xaxis = L"t", yaxis = L"\mathbf{y}(t)", label=L"\tilde{y}_%$(i)(t)",dpi=500, alpha=1.0)
    end
    
    # save the picture
    cd("plots")
    println("Saving $(plot_name).png at $(pwd())")
    savefig(plot_name*".png")
    cd("..")

end
   