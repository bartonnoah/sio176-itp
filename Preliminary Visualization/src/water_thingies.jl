using MAT, CairoMakie, GibbsSeaWater, Dates, ColorSchemes

cd(@__DIR__)

datafile = "../../data/interpolated_vars.mat"
visdir = "../vis"

data = matread(datafile)

tempdata = data["cons_temp_denoised"][:]
saldata = data["abs_sal_denoised"][:]
depthdata = data["depth"][:]

offset = (-1, 1)
tempbounds = extrema(tempdata) .+ offset
salbounds = extrema(saldata) .+ offset

salgrid = range(salbounds...; length = 500)
tempgrid = range(tempbounds...; length = 500)
densgrid = gsw_sigma0.(salgrid, tempgrid')

fig = Figure()
ax = Axis(fig[1, 1],
    title = "Water Mass Plot",
    ylabel = "Conservative Temperature (ºC)",
    xlabel = "Absolute Salinity (g/kg)"
)
scatterdata = scatter!(saldata, tempdata; color = depthdata)
contour!(salgrid, tempgrid, densgrid; color = :black)
Colorbar(fig[1, 1][1, 2], scatterdata; flipaxis = true, label = "Depth (m)")

display(fig)
save(joinpath(visdir, "density_scatterplot.png"), fig)

