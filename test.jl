using Tk
using Compat; import Compat.String

w = Toplevel()

mb = Menu(w)            ## makes menu, adds to top-level window
fmenu = menu_add(mb, "File")
omenu = menu_add(mb, "Options")

menu_add(fmenu, "Open file...", (path) -> println("Open file dialog, ..."))
menu_add(fmenu, Separator(w))   ## second argument is Tk_Separator instance
menu_add(fmenu, "Close window", (path) -> destroy(w))

cb = Checkbutton(w, "Something visible")
set_value(cb, true)     ## initialize
menu_add(omenu, cb)     ## second argument is Tk_Checkbutton instance

menu_add(omenu, Separator(w))   ## put in a separator

rb = Radio(w, ["option 1", "option 2"])
set_value(rb, "option 1")   ## initialize
menu_add(omenu, rb)     ## second argument is Tk_Radio instance

function callback(path)
  vals = map(get_value, (cb, rb))
  println(vals)
end
