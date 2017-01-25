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

F_saddle = alpha -> function (t,x)
  a = 1
  b = 1
  alpha = 1

  F1 = (x1, x2) ->
    (a * alpha) + (b*(x1^2))
  F2 = (x1, x2) ->
    (-x2)
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end

F_transcritical = alpha -> function (t,x)
  a = 1
  b = 1
  #alpha = -1

  F1 = (x1, x2) ->
    (a * alpha * x) + (b*(x1^2))
  F2 = (x1, x2) ->
    (-x2)
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end

#Supercritical pitchfork bifurcation
F_super_pitch = alpha -> function (t,x)
  a = 1
  b = -1
  #alpha = 0

  F1 = (x1, x2) ->
    (a * alpha * x1) + (b*(x1^3))
  F2 = (x1, x2) ->
    (-x2)
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end

#Subcritical pitchfork bifurcation
F_sub_pitch = alpha -> function (t,x)
  a = 1
  b = 1
  #alpha = -1

  F1 = (x1, x2) ->
    (a * alpha * x1) + (b*(x1^3))
  F2 = (x1, x2) ->
    (-x2)
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end

F_super_hopf = alpha -> function (t,x)
  #alpha = 0

  F1 = (x1, x2) ->
    -x2 + (x1 * (alpha -((x1^2)+(x2^2))))
  F2 = (x1, x2) ->
    x1 + (x2 * (alpha -((x1^2)+(x2^2))))
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end

F_sub_hopf = alpha -> function (t,x)
  #alpha = 0

  F1 = (x1, x2) ->
    -x2 + (x1 * (alpha +((x1^2)+(x2^2))))
  F2 = (x1, x2) ->
    x1 + (x2 * (alpha +((x1^2)+(x2^2))))
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end
