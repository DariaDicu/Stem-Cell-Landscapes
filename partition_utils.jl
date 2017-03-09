# You must open this URL for the code to work:
# https://www.youtube.com/watch?v=pZ12_E5R3qc

# Keep these functions in partition_utils.jl and use
# include("partition_utils.jl") to import them.
using GLVisualize, GLWindow, Colors, GeometryTypes

# Fix the right area along the x-axis to be "right_width" pixels wide.
function x_partition_fixed_right(area, right_width)
  left_area = map(area) do a
    SimpleRectangle{Int64}(a.x, a.y, a.w-right_width, a.h)
  end
  right_area = map(area) do a
    SimpleRectangle{Int64}(a.x+a.w-right_width, a.y, right_width, a.h)
  end
  left_area, right_area
end

# Fix the left area along the x-axis to be "left_width" pixels wide.
function x_partition_fixed_left(area, left_width)
  left_area = map(area) do a
    SimpleRectangle{Int64}(a.x, a.y, left_width, a.h)
  end
  right_area = map(area) do a
    SimpleRectangle{Int64}(a.x+left_width, a.y, a.w-left_width, a.h)
  end
  left_area, right_area
end

# Fix the top area along the y-axis to be "top_height" pixels high.
function y_partition_fixed_top(area, top_height)
  bottom_area = map(area) do a
    SimpleRectangle{Int64}(a.x, a.y, a.w, a.h-top_height)
  end
  top_area = map(area) do a
    SimpleRectangle{Int64}(a.x, a.y+a.h-top_height, a.w, top_height)
  end
  bottom_area, top_area
end

# Fix the bottom area along the y-axis to be "bottom_height" pixels high.
function y_partition_fixed_bottom(area, bottom_height)
  bottom_area = map(area) do a
    SimpleRectangle{Int64}(a.x, a.y, a.w, bottom_height)
  end
  top_area = map(area) do a
    SimpleRectangle{Int64}(a.x, a.y+bottom_height, a.w, a.h-bottom_height)
  end
  bottom_area, top_area
end

#=
# Example use case:
window = glscreen()

editarea, viewarea = y_partition_fixed_bottom(window.area, 180)

edit_screen = Screen(
    window, area = editarea,
    color = RGBA{Float32}(0.0f0, 0.0f0, 0.0f0, 1f0))
view_screen = Screen(
    window, area = viewarea,
    color = RGBA(255.0f0, 255.0f0, 255.0f0, 1f0),
    stroke = (1f0, RGBA{Float32}(0.13f0, 0.13f0, 0.13f0, 13f0)))

renderloop(window)
=#
