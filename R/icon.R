# - not run - #

# ## create icon
#
# img <- readImage(files = file.path(getwd(), "extras","icon-source.jpg"))
# colorMode(img) <- Grayscale
# imgt <- img < 0.9
#
# # create owin
# w <- as.polygonal(owin(mask= imageData(imgt)[,,1]))
# w <- rotate(w, -pi/2)
#
# # populate with random points
# csr <- rpoint(n = 1000, win = w)
# csr$marks <- rnorm(1000)
#
#
# svg("icon.svg", width = 7, height = 4)
# # plot.ppp(csr, show.window = FALSE, use.marks = F, cols = "cornflowerblue",
# #                    legend = F, size = 0.7, main = "", pch = 1, show.all = FALSE)
# plot.ppp(csr, show.window = FALSE, use.marks = F, cols = "black", bg = "white",
#                    legend = F, size = 0.7, main = "", pch = 21, show.all = FALSE)
# dev.off()
