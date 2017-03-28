using PyDSTool, DataStructures

name = "Calcium channel model"
pars = Dict{String,Any}(
          "vl"=> -60,
          "vca"=> 120,
          "i" => 0,
          "gl" => 2,
          "gca" => 4,
          "c" => 20,
          "v1" => -1.2,
          "v2" => 18)
vars = Dict{String,Any}(
          "v" => "( i + gl * (vl - v) - gca * 0.5 * (1 + tanh( (v-v1)/v2 )) * (v-vca) )/c",
          "w" => "v-w")
ics = Dict{String,Any}(
          "v" => 0,
          "w" => 0)
tdomain = [0;30]


dsargs = build_ode(name,ics,pars,vars,tdomain)

#Solve the ODE
#d = solve_ode(dsargs)
#using Plots
#plot(d[:t],d[:v])

#Bifurcation Plots
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
p = plot(d[:i][1:150],d[:v][1:150],color=:blue,leg=false,lw=3)
plot!(p,d[:i][150:280],d[:v][150:280],lw=3,line=(:dash),color=:red)
plot!(p,d[:i][280:end],d[:v][280:end],color=:blue,lw=3)
scatter!(p,[special_points["LP1"][:i]],[special_points["LP1"][:v]],label="LP1",markersize=15,color=:red)
scatter!(p,[special_points["LP2"][:i]],[special_points["LP2"][:v]],label="LP1",markersize=15,color=:red)
=#

d,stab,special_points = bifurcation_curve(PC,"LP-C",["i","gca"],
                              max_num_points=200,initpoint="EQ1:LP2",
                              max_stepsize=2,min_stepsize=1e-5,
                              stepsize=2e-2,loc_bif_points="CP",
                              save_eigen=true,name="SN1",
                              print_info=true,calc_stab=true,
                              solver_sequence=[:forward,:backward])

#=
using Plots
p = plot(d[:i],d[:gca],color=:blue,leg=false,lw=3)
sps = special_points["CP1"][1:2]
scatter!(p,[sps[2]],[sps[1]],color=:red,markersize=15)
=#
