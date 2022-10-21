using PyCall
using .Moa
using CSV
using LinearAlgebra: norm

function MoaPlot(state, exaggeration::Float64, toVisualize::Vector{Float64})
    # Import pyvista
    pv = PyCall.pyimport_conda("pyvista", "pyvista", "conda-forge")
    PyCall.pyimport_conda("matplotlib.pyplot", "matplotlib", "conda-forge")

    # Plot positions
    point_cloud = pv.PolyData([node.position+(node.displacement*exaggeration) for node in state.nodes])

    # Create plotter
    plotter = pv.Plotter()
    plotter.add_axes(interactive=true)
    plotter.enable()

    # Add point cloud to plotter
    plotter.add_mesh(point_cloud, scalars=toVisualize, cmap="plasma")

    # Show the plot
    plotter.show()
end

function plotDamage(state, exaggeration::Float64)
    MoaPlot(state, exaggeration, Moa.MoaUtil.GetDamageVector(state.nodes))
end

function plotInterfaceDamage(state, exaggeration::Float64)
    MoaPlot(state, exaggeration, [Float64(Moa.Nodes.interfaceDamage(node)) for node in state.nodes])
end

function plotMaterialDamage(state, exaggeration::Float64)
    MoaPlot(state, exaggeration, [Float64(Moa.Nodes.materialDamage(node)) for node in state.nodes])
end

function plotDisplacement(state, exaggeration::Float64, axis::Int64)
    MoaPlot(state, exaggeration, [node.displacement[axis] for node in state.nodes])
end

function plotDisplacementMagnitude(state, exaggeration::Float64)
    MoaPlot(state, exaggeration, [norm(node.displacement) for node in state.nodes])
end

function plotVelocity(state, exaggeration::Float64, axis::Int64)
    MoaPlot(state, exaggeration, [node.velocity[axis] for node in state.nodes])
end

function plotMaterialID(state, exaggeration::Float64)
    MoaPlot(state, exaggeration, [Float64(node.material.id) for node in state.nodes])
end

function plotVelocityMagnitude(state, exaggeration::Float64)
    MoaPlot(state, exaggeration, [norm(node.velocity) for node in state.nodes])
end


function plotOutput(filepath::String, exaggeration::Float64)
    grid = CSV.File(filepath, stripwhitespace=true,  comment="#")

    @assert :x in grid.names
    @assert :y in grid.names
    @assert :z in grid.names

    positions = fill(zeros(3), length(grid))

    for (i,row) in enumerate(grid)
        positions[i] = [Float64(row[:x]),Float64(row[:y]),Float64(row[:z])]
    end

    if :ux in grid.names
        @assert :ux in grid.names
        @assert :uy in grid.names
        @assert :uz in grid.names

        for (i,row) in enumerate(grid)
            positions[i] += [Float64(row[:ux]),Float64(row[:uy]),Float64(row[:uz])] * exaggeration
        end
    end

    # Import pyvista
    pv = PyCall.pyimport_conda("pyvista", "pyvista", "conda-forge")
    PyCall.pyimport_conda("matplotlib.pyplot", "matplotlib", "conda-forge")

    
    # Plot positions
    point_cloud = pv.PolyData(positions)

    # Create plotter
    plotter = pv.Plotter()
    plotter.add_axes(interactive=true)
    plotter.enable()

    # Add point cloud to plotter
    if :dmg in grid.names
        plotter.add_mesh(point_cloud, scalars=[Float64(row[:dmg]) for row in grid], cmap="plasma")
    else
        plotter.add_mesh(point_cloud, cmap="plasma")
    end

    # Show the plot
    plotter.show()
end
