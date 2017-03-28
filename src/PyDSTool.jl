module PyDSTool

using PyCall, DataStructures, DiffEqBase, RecipesBase

const ds = PyNULL()

function __init__()
    copy!(ds, pyimport("PyDSTool"))
end

include("constants.jl")
include("ode_construct_solve.jl")
include("bifurcation.jl")

export ds, build_args

export PYDSTOOL_CURVE_CLASSES, ALL_POINT_TYPES

export set_name,set_ics,set_pars, set_vars, set_tdata,set_fnspecs,
       set_tdomain, interpert_pts, build_ode, solve_ode,bifurcation_curve

export find_changes
end
