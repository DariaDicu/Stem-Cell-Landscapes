using Plots, ImageMagick

plt = plot([sin,cos], 0:0.1:2π)

@gif for i = 0:0.1:2π
  plot([sin,cos], 0+i:0.1:2π+i)
end every 3
