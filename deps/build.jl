using PyCall
using Conda


if PyCall.conda
	Conda.add("pip")
	pip = joinpath(Conda.BINDIR, "pip")
	run(`$pip install pydstool`)
else
	try
		pyimport("pydstool")
	catch ee
		typeof(ee) <: PyCall.PyError || rethrow(ee)
		warn("""
Python Dependancies not installed
Please either:
 - Rebuild PyCall to use Conda, by running in the julia REPL:
    - `ENV[PYTHON]=""; Pkg.build("PyCall"); Pkg.build("PyDSTool")`
 - Or install the depencences, eg by running pip
	- `pip install pydstool`
	"""
		)
	end

end
