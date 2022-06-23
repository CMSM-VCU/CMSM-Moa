using PyCall
using .Moa
using CSV

function plot(nodes::Vector{Moa.Nodes.Node}, exaggeration::Float64, toVisualize::Vector{Float64})
    # Import pyvista
    pv = PyCall.pyimport_conda("pyvista", "pyvista", "conda-forge")
    PyCall.pyimport_conda("matplotlib.pyplot", "matplotlib", "conda-forge")

    # Plot positions
    point_cloud = pv.PolyData([node.position+(node.displacement*exaggeration) for node in nodes])

    # Create plotter
    plotter = pv.Plotter()
    plotter.add_axes(interactive=true)
    plotter.enable()

    # Add point cloud to plotter
    plotter.add_mesh(point_cloud, scalars=toVisualize, cmap="plasma")

    # Show the plot
    plotter.show()
end

function plotDamage(nodes::Vector{Moa.Nodes.Node}, exaggeration::Float64)
    plot(nodes, exaggeration, Moa.MoaUtil.GetDamageVector(Moa.nodes))
end

function plotDisplacement(nodes::Vector{Moa.Nodes.Node}, exaggeration::Float64, axis::Int64)
    plot(nodes, exaggeration, [node.displacement[axis] for node in nodes])
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