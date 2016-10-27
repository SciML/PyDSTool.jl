using PyDSTool, PyCall

icdict = Dict("x"=>1,"y"=>0.4)
py_icdict = PyDict(icdict)
pardict = Dict("k"=>0.1,"m"=>0.5)
py_pardict = PyDict(pardict)
x_rhs = "y"
y_rhs = "-k*x/m"
vardict = Dict("x"=>x_rhs,"y"=>y_rhs)
py_vardict = PyDict(vardict)

DSargs = ds.args()
DSargs[:name] = "SHM"
DSargs[:ics] = py_icdict
DSargs[:pars] = py_pardict
DSargs[:tdata] = [0,20]
DSargs[:varspecs] = py_vardict
DS = ds.Generator.Euler_ODEsystem(DSargs)
