# Code for generating a Julia OpenGL model that plots the landscape for a fixed
# value of parameter a and animated moving "ant" trajectories onto the model
# from a subset of the simulations.
include("callback_registrar.jl")
include("ode_simulator.jl")
include("gui.jl")
include("render.jl")

using ODESimulator, KernelDensity, Interpolations, DataFrames, Reactive

# Code for running initial simulation.
# TODO: remove and start up with ODE GUI.
F = function (t,x)
  a = 0.3
  n = 4
  S = 0.5
  k = b = 1
  F1 = (x1, x2) ->
    (a*(x1^n)/(S^n + x1^n) + b*S^n/(S^n + x2^n) - k*x1)
  F2 = (x1, x2) ->
    (a*(x2^n)/(S^n + x2^n) + b*S^n/(S^n + x1^n) - k*x2)
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end

runs = 100 # Number of simulation runs
n = 2 # Number of dimensions

simulation_data = ODESimulator.build_landscape(runs, F, 2, (0,5))
data_signal = Signal(simulation_data)

# Function for resetting the ODE after "Run" is clicked. Runs the simulations
# again and updates the landscape.
function reset_ode_callback(ode_definition, dims, bounds, runs)
  global bool_click = true
  simulation_data = ODESimulator.build_landscape(
    runs, ode_definition, dims, bounds)
  data_signal
  println(simulation_data)
  push!(data_signal, simulation_data)
end

CallbackRegistrar.register_callback(:reset_ode, reset_ode_callback)

initialize_gui()
render(data_signal)
