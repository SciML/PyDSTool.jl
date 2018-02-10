__precompile__()

module PyDSTool

using PyCall, DataStructures, DiffEqBase, RecipesBase

const ds = PyNULL()

function __init__()
  try
    copy!(ds, pyimport_conda("PyDSTool", "pydstool", "conda-forge"))
    return
  catch err
    if err isa PyCall.PyError
      # A dirty hack to force importing PyDSTool:
      # https://github.com/JuliaDiffEq/PyDSTool.jl/issues/5
      py"""
      import scipy

      original_version = scipy.__version__
      try:
          scipy.__version__ = '0.9'
          import PyDSTool
      finally:
          scipy.__version__ = original_version
      """
    else
      rethrow()
    end
  end
  copy!(ds, pyimport_conda("PyDSTool", "pydstool", "conda-forge"))
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
