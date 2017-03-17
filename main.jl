# Code for generating a Julia OpenGL model that plots the landscape for a fixed
# value of parameter a and animated moving "ant" trajectories onto the model
# from a subset of the simulations.
include("callback_registrar.jl")
include("ode_simulator.jl")
include("gui.jl")
include("render.jl")

using ODESimulator

# Function for resetting the ODE after "Run" is clicked. Runs the simulations
# again and updates the landscape.
# TODO: Pass timespan.
function create_model_callback(ode_definition, dims, bounds, runs, time)
  ODEInput.destroy_gui_window()
  Renderer.reset_screens()
  simulation_data = ODESimulator.build_landscape_parallel(
    runs, ode_definition, dims, bounds, time)
  Renderer.rerender(simulation_data, dims, runs)
end

CallbackRegistrar.register_callback(:create_model, create_model_callback)

Renderer.init()
ODEInput.initialize_gui()
