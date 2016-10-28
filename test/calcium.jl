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

names,d,stab,special_points = bifurcation_curve(PC,"EP-C",["i"],
                          max_num_points=450,
                          max_stepsize=2,min_stepsize=1e-5,
                          stepsize=2e-2,loc_bif_points="all",
                          save_eigen=true,name="EQ1",
                          print_info=true,calc_stab=true)

function find_changes(stab)
  changes = Int[]
  for i in 2:length(stab)
    if stab[i]!= stab[i-1]
      push!(changes,i)
    end
  end
  changes
end

using Plots
p = plot(d["i"][1:148],d["v"][1:148],color=:blue,leg=false,lw=3)
plot!(p,d["i"][148:280],d["v"][148:280],lw=3,line=(:dash),color=:red)
plot!(p,d["i"][280:end],d["v"][280:end],color=:blue,lw=3)
scatter!(p,[d["i"][150]],[d["v"][150]],label="LP1",markersize=15,color=:red)
scatter!(p,[d["i"][281]],[d["v"][281]],label="LP2",markersize=15,color=:red)

names,d,stab,special_points = bifurcation_curve(PC,"LP-C",["i","gca"],
                              max_num_points=200,initpoint="EQ1:LP2",
                              max_stepsize=2,min_stepsize=1e-5,
                              stepsize=2e-2,loc_bif_points="CP",
                              save_eigen=true,name="SN1",
                              print_info=true,calc_stab=true,
                              solver_sequence=[:forward,:backward])

using Plots
p = plot(d["i"],d["gca"],color=:blue,leg=false,lw=3)
sps = special_points["CP1"][1:2]
scatter!(p,[sps[2]],[sps[1]],color=:red,markersize=15)
# Note: names[1] == "gca", names[2]=="i"

#=
#Manual
PC = ds.ContClass(ode)
PCargs = ds.args(name="EQ1")
PCargs[:type]         = "EP-C"
PCargs[:freepars]     = ["i"]
PCargs[:MaxNumPoints] = 450
PCargs[:MaxStepSize]  = 2
PCargs[:MinStepSize]  = 1e-5
PCargs[:StepSize]     = 2e-2
PCargs[:LocBifPoints] = "all"
#PCargs[:LocBifPoints] = "LP"
PCargs[:SaveEigen]    = true

PC[:newCurve](PCargs)
PC[:curves]["EQ1"][:forward]()

#More data
#PC[:curves]["EQ1"][:sol][148][:labels]["EP"]["data"][:__dict__]

PC[:curves]["EQ1"][:info]()

pts = PC[:curves]["EQ1"][:_curveToPointset]()
d = OrderedDict{String,Vector{Float64}}()
names = pts[:_ix_name_map]
arrs = pts[:coordarray]
for i in 1:length(names)
  d[names[i]] = arrs[i,:]
end
len = length(d[first(keys(d))])

# S => Stable
# U => Unstable
# N => Neutral
stab = [PC[:curves]["EQ1"][:sol][i][:labels]["EP"]["stab"] for i in 1:len]

# Get information for type of points, ex limit points
special_points = Dict{String,Any}()
for point_type in ALL_POINT_TYPES
  res = PC[:curves]["EQ1"][:sol][:bylabel](point_type)
  if res != nothing
    for i in 1:length(res)
      special_points[point_type*string(i)] = res[i][:coordarray]
    end
  end
end

#Return names,d,stab,special_points
=#
