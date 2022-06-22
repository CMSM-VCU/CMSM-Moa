using PyCall
using .Moa


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