using PyDSTool, PyCall

icdict = Dict("x"=>1,"y"=>0.4)
py_icdict = PyDict(icdict)
pardict = Dict("k"=>0.1,"m"=>0.5)
py_pardict = PyDict(pardict)
x_rhs = "y"
y_rhs = "-k*x/m"
vardict = Dict("x"=>x_rhs,"y"=>y_rhs)
py_vardict = PyDict(vardict)

# keys()...

DSargs = ds.args()
DSargs[:name] = "SHM"
DSargs[:ics] = py_icdict
DSargs[:pars] = py_pardict
DSargs[:tdata] = [0,20]
DSargs[:varspecs] = py_vardict
DS = ds.Generator[:Vode_ODEsystem](DSargs)
traj = DS[:compute]("demo")
pts = traj[:sample]()
d = Dict{Symbol,Vector{Float64}}()
d[Symbol(pts[:indepvarname])] = pts[:indepvararray]
names = pts[:_ix_name_map]
arrs = pts[:coordarray]
for i in 1:length(names)
  d[Symbol(names[i])] = arrs[i,:]
end
using Plots
plot(d[:t],[d[:x],d[:y]])

# Interpolations
t = 5.4
traj(t)[:coordarray]
