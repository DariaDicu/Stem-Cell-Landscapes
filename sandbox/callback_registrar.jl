# Module for registering callbacks between communicating scripts. This keeps
# the dependencies between scripts less tangled, since a script can register a
# callback and another can execute it using the name Symbol (e.g. :gui_popup).
module CallbackRegistrar
  # e.g. Dict[:gui_popup => function my_popup_code() ... end, ...]
  callbacks = Dict()
  function register_callback(name::Symbol, callback)
    global callbacks
    callbacks[name] = callback
  end

  function get_callback(name::Symbol)
    global callbacks
    callbacks[name]
  end
end
