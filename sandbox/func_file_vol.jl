(a, n, b, S, k, vol) -> begin println("Function evaluated!"); function (t,x)
    F1 = (x1, x2) ->
      (a*(x1^n)/(S^n + x1^n) + b*S^n/(S^n + x2^n) - k*x1 + vol*randn())
    F2 = (x1, x2) ->
      (a*(x2^n)/(S^n + x2^n) + b*S^n/(S^n + x1^n) - k*x2 + vol*randn())
    return [F1(x[1], x[2]), F2(x[1], x[2])]
  end
end
