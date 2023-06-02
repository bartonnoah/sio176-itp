using MAT, GibbsSeaWater, Dates, ColorSchemes, PyPlot, StatsBase

cd(@__DIR__)

datadir = "../../data"
visdir = "../vis"

data = matread(joinpath(datadir, "itp121_data.mat"))
derived_vars = matread(joinpath(datadir, "derived_variables.mat"))

milliseconds_in_day = Dates.value(Millisecond(Day(1)))

times = @. Millisecond(round(Int, data["time"] * milliseconds_in_day)) + DateTime(0,1,1)

depths = derived_vars["depth"]

cons_temp = derived_vars["cons_temp"]

pot_dens_anom = derived_vars["pot_dens_anom"]

#plots demands a square times array
square_times = repeat(times, size(depths, 1))

num_times = Float64.(Dates.value.(square_times))

skipna(x) = filter(!isnan, x)

#Now remove depth Nan

cons_temp[depths .=== NaN] .= NaN

avg_depths = Dict(axes(depths, 1).=>mapslices(x->mean(filter(!isnan, x)), depths; dims=2))

for idx in CartesianIndices(depths)
    if isnan(depths[idx])
        first_idx = first(Tuple(idx))
        depths[idx] = avg_depths[first_idx]
    end
end

#Now plot
p = contourf(num_times, depths, cons_temp, cmap = :inferno)
plt.colorbar(mappable = p, label = "Conservative Temperature (ºC)")
gca().invert_yaxis()
plt.title("Conservative Temperature with Density Overlay")
plt.xlabel("Time")
plt.ylabel("Depth(m)")

#Handle time ticks
t_ticks = range(extrema(times)...; step=(maximum(times) - minimum(times))/8)
t_tickvals = Dates.value.(t_ticks)
t_ticklabels = Dates.format.(t_ticks, dateformat"yyyy-mm-dd")
plt.xticks(ticks=t_tickvals, labels = t_ticklabels, rotation=20)

#Now contour stuff
cs = contour(num_times, depths, pot_dens_anom; colors = "grey")
label_yvals = [30, 40, 50, 120, 210]
label_pos = tuple.(Dates.value(Date(2021, 11, 28)), label_yvals)
gca().clabel(cs, cs.levels, inline=true)

plt.tight_layout()
savefig(joinpath(visdir, "pot_dens_on_cons_temp.png"))
plt.close()


