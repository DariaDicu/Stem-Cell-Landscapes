my_func = function ()
  println("Hello")
end

function executing_function(x)
  x()
end

executing_function(function ()
  println("Hello")
end)
