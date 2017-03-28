using ParameterizedFunctions, PyDSTool, Plots

f = @ode_def_bare Calcium begin
  dv = ( i + gl * (vl - v) - gca * 0.5 * (1 + tanh( (v-v1)/v2 )) * (v-vca) )/c
  dw = v-w
end vl=>-60 vca=>120 i=>0.0 gl=>2 gca=>4 c=>20 v1=>-1.2 v2=>18
u0 = [0;0]
tspan = [0;30]

dsargs = build_ode(f,u0,tspan)
#Solve the ODE
d = solve_ode(dsargs)
#using Plots
#plot(d[:t],d[:v])

ode = ds[:Generator][:Vode_ODEsystem](dsargs)
ode[:set](pars = Dict("i"=>-220))
ode[:set](ics  = Dict("v"=>-170))
PC = ds[:ContClass](ode)

d,stab,special_points = bifurcation_curve(PC,"EP-C",["i"],
                          max_num_points=450,
                          max_stepsize=2,min_stepsize=1e-5,
                          stepsize=2e-2,loc_bif_points="all",
                          save_eigen=true,name="EQ1",
                          print_info=true,calc_stab=true)

#=
using Plots, PlotRecipes
y = collect((1:200)/11)
x = -cos(y)
g = ["N"; [x[i] > x[i-1] ? "N" : "S" for i in 2:length(x)]]
compositeline(x, y, g, linestyle = [:solid :dash], color = [:blue :red])
=#
