using MAT, Plots, GibbsSeaWater, Dates, ColorSchemes

pyplot()
cd(@__DIR__)

datafile = "../../data/itp121_data.mat"
visdir = "../vis"

data = matread(datafile)

milliseconds_in_day = Dates.value(Millisecond(Day(1)))

times = @. Millisecond(round(Int, data["time"] * milliseconds_in_day)) + DateTime(0,1,1)

pressure = data["pres"]

depth = -1 .* gsw_z_from_p.(pressure, data["latitude"], 0, 0)

times_expanded = repeat(times, size(pressure, 1))

abs_sal = gsw_sa_from_sp.(data["salinity"], pressure, data["longitude"], data["latitude"])

abs_sal[abs.(abs_sal).>2000] .= NaN

cons_temp = gsw_ct_from_t.(abs_sal, data["temperature"], pressure)

pot_dens_anom = gsw_sigma0.(abs_sal, cons_temp)

varnames = ["Absolute salinity", "Conservative temperature", "Potential density anomaly"]
vardata = [abs_sal, cons_temp, pot_dens_anom]
varunits = ["g/kg", "ºC", "kg/m^3"]
colors = [:haline, :thermal, :dense]

for (name, data, unit, color) in zip(varnames, vardata, varunits, colors)
    p = contourf(times_expanded, depth, data; title = "$name vs depth and time", xlabel = "time", ylabel = "Depth (m)", colorbar_title = "$name ($unit)", yflip=true, c = color, dpi = 1200)
    display(p)
    savefig(p,joinpath(visdir, "$name.png"))
end