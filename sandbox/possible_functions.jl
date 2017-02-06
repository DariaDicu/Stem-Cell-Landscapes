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

#Quad-stable system from Zhou et al, 2012
F_quad_stable = function (t,x)

  F1 = (x1, x2) ->
   -1 + 9*x1 -2(x1^3) + 9*x2 -2(x2^3)
  F2 = (x1, x2) ->
   1 - 11*x1 + 2(x1^3) + 11*x2 -2(x2^3)
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end

#Stumpf's proposed toy system
# A ---> A  , with rate a
# A ---> B  , with rate b
# B ---| A  , with rate g
F_toy_system = function (t,x)
  a = 1
  b = 1
  g = 1

  F1 = (x1, x2) ->
   (a * x1) - (b * x1) - (g * x2)
  F2 = (x1, x2) ->
   b*x1
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end

#Simple bistability with two mutual repressors
# A ---| B  , with rate alpha * i
# B ---| A  , with rate alpha * j
F_mutual_repressors = function (t,x)
  alpha = 1
  i = 1
  j = 1

  F1 = (x1, x2) ->
   -alpha * j * x2
  F2 = (x1, x2) ->
   -alpha * i * x1
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end

#More involved bistable system
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3213676/#B13
F_official_bistable = function (t,x)

  r = 2 #Rate constants
  b = 0.2 #Effective affinity constants

  kyx = 0.7 #Effective affinity
  kxy = 0.5 #Effective affinity

  deg = 1 #First order degradation
  n = 4   #Hill coefficient

  F1 = (x1, x2) ->
   b - deg*x1 + (r*(kyx^n))/((kyx^n)+(x2^n))
  F2 = (x1, x2) ->
   b - deg*x2 + (r*(kxy^n))/((kxy^n)+(x1^n))
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end

#Saddle node bifurcation
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

#Transcitical bifurcation
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

#Supercritical Hopf bifurcation
F_super_hopf = alpha -> function (t,x)
  #alpha = 0

  F1 = (x1, x2) ->
    -x2 + (x1 * (alpha -((x1^2)+(x2^2))))
  F2 = (x1, x2) ->
    x1 + (x2 * (alpha -((x1^2)+(x2^2))))
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end

#Subcritical Hopf bifurcation
F_sub_hopf = alpha -> function (t,x)
  #alpha = 0

  F1 = (x1, x2) ->
    -x2 + (x1 * (alpha +((x1^2)+(x2^2))))
  F2 = (x1, x2) ->
    x1 + (x2 * (alpha +((x1^2)+(x2^2))))
  return [F1(x[1], x[2]), F2(x[1], x[2])]
end
