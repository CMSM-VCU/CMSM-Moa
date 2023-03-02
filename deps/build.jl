import Pkg
Pkg.add("CSV")
Pkg.add("StaticArrays")
Pkg.add("Revise")

# PyCall
Pkg.add("PyCall")
# Pycall does not use julia installation by default in linux
if Sys.islinux()
    ENV["PYTHON"] = ""
    Pkg.build("PyCall")
end

using PyCall
PyCall.pyimport_conda("tables", "pytables")
PyCall.pyimport_conda("pandas", "pandas")

# Vislualization
# (not for headless computers)
if Sys.iswindows()
    PyCall.pyimport_conda("pyvista", "pyvista", "conda-forge")
    PyCall.pyimport_conda("matplotlib.pyplot", "matplotlib", "conda-forge")
    Pkg.add("Plots")
end
