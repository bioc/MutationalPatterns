default_colors_ggplot <- function(n) {
  
  hues = seq(15, 375, length = n + 1)
  return(hcl(h = hues, l = 65, c = 100)[1:n])

}

