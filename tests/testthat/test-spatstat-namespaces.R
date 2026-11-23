# These tests run without attaching any spatstat package. They exercise the
# geometry, image, simulation, and plotting paths formerly supplied by import.

test_that("the package namespace contains no dynamic spatstat copies", {
    namespace <- asNamespace("moleculaR")

    for (name in c("installedPckgs", "spatstatFuncs", "correspondPckgs", "ppp", "density.ppp",
        "summary.owin", "as.polygonal")) {
        expect_false(exists(name, envir = namespace, inherits = FALSE))
    }

    expect_false("package:spatstat" %in% search())
    expect_true(is.function(spatstat.geom::summary.owin))
})

test_that("spatial combination, image conversion, and expression evaluation work", {
    spp <- searchAnalyte(100, 0.1, fixture_cardinal())
    image <- spp2im(spp, rescale = FALSE)
    roundtrip <- im2spp(image, win = spp$window, rescale = FALSE)

    expect_s3_class(image, "im")
    expect_s3_class(roundtrip, "analytePointPattern")
    expect_equal(sort(roundtrip$marks$intensity), sort(spp$marks$intensity))

    combined <- superimposeAnalytes(list(spp, spp), spWin = spp$window, check = FALSE)
    expect_equal(combined$n, 2 * spp$n)

    selected <- subsetAnalytes(combined, mzVals == 100)
    expect_equal(selected$n, combined$n)

    expression <- evalSpatialExpr("a + a", dataList = list(a = spp), win = spp$window)
    expect_s3_class(expression, "analytePointPattern")
    expect_equal(sort(expression$marks$intensity), sort(2 * spp$marks$intensity))

    # Two squares of area four overlap by area two: DSC = 2*2 / (4+4).
    first <- spatstat.geom::owin(c(0, 2), c(0, 2))
    second <- spatstat.geom::owin(c(1, 3), c(0, 2))
    expect_equal(dsc(first, first), 1)
    expect_equal(dsc(first, second), 0.5)

    expect_true(all(is.finite(.bw.scott2(spp))))
    expect_true(is.finite(.bw.scott.iso2(spp)))
})

test_that("all simulation shapes resolve their geometry explicitly", {
    circle <- createHotspotWin("centralCircle", radius = 2, backgroundLength = 20)
    multiple <- createHotspotWin("multiCircle", radius = 2, backgroundLength = 20)
    donut <- createHotspotWin("donut", radius = c(1, 2), backgroundLength = 20)

    expect_s3_class(circle, "owin")
    expect_s3_class(multiple, "owin")
    expect_s3_class(donut, "owin")
    expect_gt(spatstat.geom::area.owin(multiple), spatstat.geom::area.owin(circle))
    expect_lt(spatstat.geom::area.owin(donut), spatstat.geom::area.owin(circle))

    simulated <- simulatespp(lambdaHotspot = 1, lambdaBackground = 0.1, intensityData = seq_len(10000),
        radius = 2, backgroundLength = 20, seed = 42, sppPlot = FALSE)

    expect_s3_class(simulated, "analytePointPattern")
    expect_gt(simulated$n, 0)
})

test_that("plotting and smoothing do not rely on attached spatial packages", {
    spp <- searchAnalyte(100, 0.1, fixture_cardinal())
    mpm <- suppressWarnings(probMap(spp, bw = 1, seed = 1, verbose = FALSE))

    # A file device makes these smoke tests usable on headless CI machines.
    grDevices::pdf(tempfile(fileext = ".pdf"), width = 8, height = 8)
    on.exit(grDevices::dev.off(), add = TRUE)

    expect_no_error(plotAnalyte(spp))
    expect_no_error(plotImg(spp2im(spp), smooth = TRUE))
    expect_no_error(plot(mpm, what = "MPM", smooth = TRUE, densityAsBase = FALSE, figLegend = FALSE))
})
