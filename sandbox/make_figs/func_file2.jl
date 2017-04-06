# Example input function for gr_gif_generator.jl.
(t,x) -> begin
  N,O,F,G = x
  k0=0.005
  c0=0.01
  c1=0.4
  c2=1.0
  c3=0.1
  c4=0.00135
  v=0.01
  a=0.95
  e0=0.01
  e1=1.0
  e2=1.0
  a0=0.01
  a1=1.0
  a2=5.0
  b0=0.005
  b1=0.005
  b2=1.0
  b3=1.0
  LL=0.0
  II=0.0
  dN=k0*O*(c0+c1*N*N+k0*O+c2*LL)/(1+k0*O*(c0+c1*N^2+k0*O+c2*LL+c3*F*F)+c4*O*G*G)-v*N
  dO=a+(e0+e1*O)/(1.0+e1*O+e2*G*G)-v*O
  dF=(a0+a1*O)/(1.0+a1*O+a2*II)-v*F
  dG=(b0+b1*G*G+b3*O)/(1.0+b1*G*G+b2*N*N+b3*O)-v*G
  [dN, dO, dF, dG]
end
