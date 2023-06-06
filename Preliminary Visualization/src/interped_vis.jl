using MAT, Plots, GibbsSeaWater, Dates, ColorSchemes

pyplot()
cd(@__DIR__)

datafile = "../../data/interpolated_vars.mat"
visdir = "../vis"

data = matread(datafile)

milliseconds_in_day = Dates.value(Millisecond(Day(1)))

times = @. Millisecond(round(Int, data["time"] * milliseconds_in_day)) + DateTime(0,1,1)

varnames = ["absolute salinity", "conservative temperature", "potential density anomaly"]
vardata = [data["abs_sal_denoised"], data["cons_temp_denoised"], data["pot_dens_anom_denoised"]]
varunits = ["g/kg", "ºC", "kg/m^3"]
colors = [:haline, :thermal, :dense]

depth=data["depth"]
times_expanded = repeat(times, size(depth, 1))

for (name, data, unit, color) in zip(varnames, vardata, varunits, colors)
    p = contourf(times_expanded, depth, data; title = "Smoothed $name vs depth and time", xlabel = "time", ylabel = "Depth (m)", colorbar_title = "$name ($unit)", yflip=true, c = color, dpi = 1200)
    display(p)
    savefig(p,joinpath(visdir, "interp_$name.png"))
end