import Pkg
Pkg.add("CSV")
Pkg.add("StaticArrays")


# If on linux, you need to run:
    # ENV["PYTHON"] = ""
    # Pkg.build("PyCall")
# once so that it will be able to automatically install packages

Pkg.add("PyCall")
using PyCall
PyCall.pyimport_conda("tables", "tables")
PyCall.pyimport_conda("pandas", "pandas")


# Dont add if using headless:
if Sys.iswindows()
    PyCall.pyimport_conda("pyvista", "pyvista", "conda-forge")
    PyCall.pyimport_conda("matplotlib.pyplot", "matplotlib", "conda-forge")
    Pkg.add("Plots")
end
