# Example input function for gr_gif_generator.jl.
(a, n, S, k, b) -> begin println("Function evaluated!"); function (t,x)
    #a = 0.3
    #n = 4
    #S = 0.5
    #k = b = 1
    F1 = (x1, x2) ->
      (a*(x1^n)/(S^n + x1^n) + b*S^n/(S^n + x2^n) - k*x1)
    F2 = (x1, x2) ->
      (a*(x2^n)/(S^n + x2^n) + b*S^n/(S^n + x1^n) - k*x2)
    #println("Function evaluated inside!")
    return [F1(x[1], x[2]), F2(x[1], x[2])]
  end
end
