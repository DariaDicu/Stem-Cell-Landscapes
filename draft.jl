my_dict = Dict()

my_dict[:hi] = "hello"


nothin = function() println("hello") end

x=Dict()

x[:hi] = nothin

typeof(x)
x[:hi] = nothin
x[:hi]

module MyModule2
  x = 112
  function reset_x(_x)
    global x = _x
  end
end


modul
