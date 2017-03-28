module PyDSTool

using PyCall, DataStructures, DiffEqBase

const ds = PyNULL()

function __init__()
    copy!(ds, pyimport("PyDSTool"))
end

PYDSTOOL_CURVE_CLASSES = Set(["EP-C","LP-C","H-C1","H-C2","FP-C","LC-C"])

CONT_BIF_POINTS = ["B", "SP"]
EQUILIBRIUM_BIF_POINTS = ["BP", "LP", "H"]
FOLD_BIF_POINTS = ["BT", "ZH", "CP"]
HOPF_BIF_POINTS = ["BT", "ZH", "GH", "DH"]
FIXEDPOINT_BIF_POINTS = ["BP", "PD", "LPC", "NS"]
FOLD_MAP_BIF_POINTS = ["CP"]
LIMITCYCLE_BIF_POINTS = ["PD", "LPC", "NS"]
OTHER_SPECIAL_POINTS = ["RG", "UZ", "P", "MX", "B"]

ALL_POINT_TYPES = union([CONT_BIF_POINTS;EQUILIBRIUM_BIF_POINTS;FOLD_BIF_POINTS;
                   HOPF_BIF_POINTS;FIXEDPOINT_BIF_POINTS;FOLD_MAP_BIF_POINTS;
                   LIMITCYCLE_BIF_POINTS;OTHER_SPECIAL_POINTS])

set_name(dsargs,name::String) = (dsargs[:name] = name; nothing)

set_ics(dsargs,icdict::PyDict) = (dsargs[:ics] = icdict; nothing)
set_ics{T}(dsargs,icdict::Dict{String,T}) = set_ics(dsargs,PyDict(icdict))

set_pars(dsargs,pardict::PyDict) = (dsargs[:pars] = pardict; nothing)
set_pars{T}(dsargs,pardict::Dict{String,T}) = set_pars(dsargs,PyDict(pardict))

set_vars(dsargs,vardict::PyDict) = (dsargs[:varspecs] = vardict; nothing)
set_vars{T}(dsargs,vardict::Dict{String,T}) = set_vars(dsargs,PyDict(vardict))

set_fnspecs(dsargs,specsdict::PyDict) = (dsargs[:fnspecs] = specsdict; nothing)
set_fnspecs{T}(dsargs,specsdict::Dict{String,T}) = set_vars(dsargs,PyDict(specsdict))

set_tdata(dsargs,tdata) = (dsargs[:tdata] = tdata; nothing)
set_tdomain(dsargs,tdomain) = (dsargs[:tdomain] = tdomain; nothing)

function build_ode(name,ics,pars,vars,tdomain)
  dsargs = ds[:args]()
  set_name(dsargs,name)
  set_ics(dsargs,ics)
  set_pars(dsargs,pars)
  set_vars(dsargs,vars)
  set_tdomain(dsargs,tdomain)
  dsargs
end

function build_ode(f::AbstractParameterizedFunction,u0,tspan)
  name = string(typeof(f))
  pars = Dict{String,Any}(); vars = Dict{String,Any}(); ics = Dict{String,Any}()
  for i in 1:length(f.params)
    pars[string(f.params[i])] = getfield(f,f.params[i])
  end
  for i in 1:length(f.syms)
    vars[string(f.syms[i])] = string(f.funcs[i])
    ics[string(f.syms[i])] = u0[i]
  end
  build_ode(name,ics,pars,vars,tspan)
end

function solve_ode(dsargs,alg=:Vode_ODEsystem,name="Default Name")
  DS = ds[:Generator][alg](dsargs)
  traj = DS[:compute](name)
  d = interpert_traj(traj)

  #=
  # Interpolations
  t = 5.4
  traj(t)[:coordarray]
  =#
end

function interpert_traj(traj)
  d = Dict{Symbol,Vector{Float64}}()
  d[Symbol(traj[:indepvarname])] = first(first(values(traj[:underlyingMesh]())))
  depvars = PyDict(traj[:sample]())
  for k in keys(depvars)
    d[Symbol(k)] = depvars[k]
  end
  d
end

function bifurcation_curve(PC,bif_type,freepars;max_num_points=450,
                          max_stepsize=2,min_stepsize=1e-5,
                          stepsize=2e-2,loc_bif_points="all",
                          save_eigen=true,name="DefaultName",
                          print_info=true,calc_stab=true,
                          initpoint=nothing,solver_sequence=[:forward])

  curve_point_type = bif_type[1:end-2]

  if !(typeof(freepars)<:AbstractArray)
    freepars = [freepars]
  end

  # Setup Parameters
  PCargs = ds[:args](name=name)
  PCargs[:type]         = bif_type
  PCargs[:freepars]     = freepars
  PCargs[:MaxNumPoints] = max_num_points
  PCargs[:MaxStepSize]  = max_stepsize
  PCargs[:MinStepSize]  = min_stepsize
  PCargs[:StepSize]     = stepsize
  PCargs[:LocBifPoints] = loc_bif_points
  PCargs[:SaveEigen]    = save_eigen
  if initpoint != nothing
    PCargs[:initpoint] = initpoint
  end

  # Run Solver
  PC[:newCurve](PCargs)
  for step in solver_sequence
    PC[:curves][name][step]()
  end

  # Print Info
  if print_info
    PC[:curves][name][:info]()
  end


  # Get the curve
  pts = PyDict(PC[:curves][name][:_curveToPointset]())
  d = OrderedDict{Symbol,Vector{Float64}}()
  for k in keys(pts)
    d[Symbol(k)] = pts[k]
  end
  len = length(d[first(keys(d))])

  # Get the stability
  # S => Stable
  # U => Unstable
  # N => Neutral

  # Get this from the information at:
  # https://github.com/robclewley/pydstool/blob/master/PyDSTool/PyCont/ContClass.py#L218
  curve = PC[:curves][name]
  if calc_stab
    stab = [curve[:CurveInfo][i][curve_point_type]["stab"] for i in 1:len]
  else
    stab = []
  end

  # Get information for special points, ex limit points
  special_points = Dict{String,Any}()
  for k in keys(curve[:BifPoints])
    for i in 1:length(curve[:BifPoints][k][:found])
      tmp_dict = PyDict(curve[:BifPoints][k][:found][i]["X"])
      dd = Dict{Symbol,Float64}()
      for k2 in keys(tmp_dict)
        dd[Symbol(k2)]=tmp_dict[k2]
      end
      special_points[k*string(i)] = dd
    end
  end
  #=
  # Start and endpoints
  for i in 1:len
    if "P" in keys(curve[i])
      tmp_dict = PyDict(curve[i]["P"]["data"]["V"])
      dd = Dict{Symbol,Float64}()
      for k in keys(tmp_dict)
        dd[Symbol(k)]=tmp_dict[k]
      end
      special_points[curve[i]["P"]["name"]] = dd
    end
  end
  =#
  d,stab,special_points
end

function find_changes(stab)
  changes = Int[]
  for i in 2:length(stab)
    if stab[i]!= stab[i-1]
      push!(changes,i)
    end
  end
  changes
end

export ds, build_args

export PYDSTOOL_CURVE_CLASSES, ALL_POINT_TYPES

export set_name,set_ics,set_pars, set_vars, set_tdata,set_fnspecs,
       set_tdomain, interpert_pts, build_ode, solve_ode,bifurcation_curve

export find_changes
end
