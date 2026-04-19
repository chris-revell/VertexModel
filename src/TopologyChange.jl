#
#  TopologyChange.jl
#  VertexModel
#
#  Created by Christopher Revell on 01/02/2021.
#
#
# Function to recalculate derived matrices (Ā, Aᵀ etc.) for any change in vertex network topology (ie any change to A or B matrices).

module TopologyChange

# Julia packages
using LinearAlgebra
using UnPack
using SparseArrays
using FastBroadcast
using FromFile
using DrWatson

# Local modules
@from "SenseCheck.jl" using SenseCheck
@from "OrderAroundCell.jl" using OrderAroundCell

function topologyChange!(R,params,matrices)

    @unpack initialSystem, 
        boundaryType,
        nCells,
        nVerts,
        nEdges,
        Λ_AA,
        Λ_AB,
        Λ_BB,
        Λ_AE,
        Λ_BE = params

    @unpack A,
        B,
        Aᵀ,
        Ā,
        Āᵀ,
        Bᵀ,
        B̄,
        B̄ᵀ,
        C,
        cellEdgeCount,
        cellVertexOrders,
        cellEdgeOrders,
        boundaryVertices,
        boundaryEdges,
        boundaryCells,
        cellPositions,
        cellLabels,
        Λs = matrices

    # Find adjacency matrices from incidence matrices
    @.. thread = false Ā .= abs.(A)    # All -1 components converted to +1 (In other words, create adjacency matrix Ā from incidence matrix A)
    @.. thread = false B̄ .= abs.(B)    # All -1 components converted to +1 (In other words, create adjacency matrix B̄ from incidence matrix B)

    # C adjacency matrix. Rows => cells; Columns => vertices. C .= B̄*Ā.÷2 (NB Integer division)
    C .= B̄ * Ā ./2

    # Update transpose matrices
    Aᵀ .= sparse(transpose(A))
    Āᵀ .= abs.(Aᵀ)
    Bᵀ .= sparse(transpose(B))
    B̄ᵀ .= abs.(Bᵀ)

    dropzeros!(A)
    dropzeros!(B)
    dropzeros!(C)
    dropzeros!(Ā)
    dropzeros!(B̄)
    dropzeros!(Aᵀ)
    dropzeros!(Āᵀ)
    dropzeros!(Bᵀ)
    dropzeros!(B̄ᵀ)

    # Calculate additional topology data
    # Number of edges around each cell found by summing columns of B̄
    cellEdgeCount .= sum.(eachrow(B̄))  # FastBroadcast doesn't work for this line; not sure why 

    # Only do the following if the initialSystem isn't periodic: 
    if boundaryType == "free"
        # Find boundary vertices
        # Summing each column of B finds boundary edges (for all other edges, cell orientations on either side cancel);
        # multiplying by Aᵀ gives nonzero values only where a vertex (row) has nonzero values at columns (edges) corresponding to nonzero values in the list of boundary edges.
        # Note that the abs is needed in case the direction of boundary edges cancel at a vertex
        boundaryVertices .= Āᵀ * abs.(sum.(eachcol(B))) .÷ 2

        # Find list of edges at system periphery
        boundaryEdges .= abs.([sum(x) for x in eachcol(B)])
    end

    # Define vectos of Λs:
    for j in 1:nEdges

        sj = dot(matrices.B̄[:, j], cellLabels)  # sparse matrix multiplication
        # In the free boundary case, check whether edge is on the boundary:
        if boundaryType=="free" && boundaryEdges[j] == 1
            if sj == 0
                Λs[j] = Λ_AE       # edge between A cell and external environment
            elseif sj == 1
                Λs[j] = Λ_BE       # edge between B cell and external environment
            else
                error("Error: edge with more than one cell on either side")
            end
        else

            if sj == 0
                Λs[j] = Λ_AA       # edge between two A cells
            elseif sj == 1
                Λs[j] = Λ_AB       # edge between A and B cell
            else
                Λs[j] = Λ_BB       # edge between two B cells
            end
        end

        
    end


    for i = 1:length(cellVertexOrders)
        cellVertexOrders[i], cellEdgeOrders[i] = orderAroundCell(matrices, i)
    end

    return nothing

end

export topologyChange!

end
