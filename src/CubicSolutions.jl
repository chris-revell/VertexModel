# Code to plot the solutions to the energy minimisation cubic used to select the number of 
# cells for the domain


module CubicSolutions
using Polynomials
using Roots
using Plots

function cubic_coeffs(Γ, Λ)
    # Function to define the cubic coefficients for an isolated hexagon of side l, for energy minimisation
    a = 9/2
    b = 0.0
    c = 12*Γ - 3*sqrt(3)
    d = Λ
    return a, b, c, d
end

function cubic_roots(Γ, Λ)
    # Function to compute the roots of the cubic for given Γ and Λ
    a, b, c, d = cubic_coeffs(Γ, Λ)
    p = Polynomial([d, c, b, a])  # Coefficients are in ascending order
    roots = roots(p)
    return roots

end

function num_real_roots(Γ,Λ)
    # Count the number of real roots of the cubic 
    roots = cubic_roots(Γ, Λ)
    return count(isreal, roots)
end

function plot_parameter_space(len)
    Λs = range(-2.0,0.0,length = len)
    Γs = range(0.0,2.0,length = len)

    N = [num_real_roots(Γ, Λ) for Λ in Λs, Γ in Γs]

    plt = heatmap(
    Γs,
    Λs,
    N,
    xlabel="Γ",
    ylabel="Λ",
    title="Number of Real Roots",
    colorbar_ticks=([1, 2, 3], ["1", "2", "3"]),
    aspect_ratio=:equal
    )

    display(plt)
    return plt
end

export plot_parameter_space

    
end
