# Precompiled vignette to speed up installation of package.
# Must manually move image files from "moleculaR/figure" direcotry to "moleculaR/vignettes/figure" after knit

knitr::knit(input = file.path(getwd(), "vignettes/moleculaR-walkthrough.Rmd.orig"),
            output = file.path(getwd(), "vignettes/moleculaR-walkthrough.Rmd"))

# more info:
# https://blog.r-hub.io/2020/06/03/vignettes/
# https://ropensci.org/blog/2019/12/08/precompute-vignettes/

