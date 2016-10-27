using PyCall
using Conda


if PyCall.conda
	Conda.add("pip")
	if is_windows() # Windows needs scipy and uses the wrong pip location
    Conda.add("scipy")
    run(`$(joinpath(Conda.PYTHONDIR, "python")) -m pip install pydstool`)
  else
    pip = joinpath(Conda.BINDIR, "pip")
    run(`$pip install pydstool`)
  end 
else
	try
		pyimport("pydstool")
	catch ee
		typeof(ee) <: PyCall.PyError || rethrow(ee)
		warn("""
Python Dependancies not installed
Please either:
 - Rebuild PyCall to use Conda, by running in the julia REPL:
    - `ENV["PYTHON"]=""; Pkg.build("PyCall"); Pkg.build("PyDSTool")`
 - Or install the depencences, eg by running pip
	- `pip install pydstool`
	"""
		)
	end

end
