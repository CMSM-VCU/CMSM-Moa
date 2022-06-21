using PyCall
using .Moa.Nodes

function plot(nodes::Vector{Moa.Nodes.Node})
    # Plotting with pyvista
    pv = PyCall.pyimport_conda("pyvista", "pyvista", "conda-forge")
    PyCall.pyimport_conda("matplotlib.pyplot", "matplotlib", "conda-forge")
    point_cloud = pv.PolyData([node.position+(node.displacement*100) for node in nodes])
    plotter = pv.Plotter()
    plotter.add_mesh(point_cloud, scalars=[node.displacement[1] for node in nodes], cmap="plasma")
    plotter.show()
end
