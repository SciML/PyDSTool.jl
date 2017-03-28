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
