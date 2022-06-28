import Pkg
Pkg.add("CSV")
Pkg.add("StaticArrays")

Pkg.add("PyCall")
PyCall.pyimport_conda("pytables", "pytables")
PyCall.pyimport_conda("pyvista", "pyvista", "conda-forge")
PyCall.pyimport_conda("matplotlib.pyplot", "matplotlib", "conda-forge")

# Dont add if using headless:
# Pkg.add("Plots")
