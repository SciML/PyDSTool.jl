using PyCall
using Conda


if PyCall.conda
    Conda.add("scipy")
    Conda.add_channel("conda-forge")
    Conda.add("pydstool")
else
	try
		pyimport("PyDSTool")
	catch ee
		typeof(ee) <: PyCall.PyError || rethrow(ee)
		warn("""
				Python Dependancies not installed
				Please rebuild PyCall to use Conda, by running in the julia REPL:
				    - `ENV["PYTHON"]=""; Pkg.build("PyCall"); Pkg.build("PyDSTool")`
				 """)
	end

end
