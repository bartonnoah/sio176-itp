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

abs_sal = derived_vars["abs_sal"]

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


diff_array = copy(pot_dens_anom)
mixed_depths = zeros(size(diff_array, 2))

for (i,col) in enumerate(eachcol(diff_array))
    col .-= first(skipna(col))
    mixed_idx = findfirst(>(0.1), col)
    if (mixed_idx) === nothing || (count(isnan, col[begin:mixed_idx]) >= 20)
        mixed_depths[i] = NaN
        continue
    end
    mixed_depths[i] = depths[mixed_idx, i]
end

p = contourf(num_times, depths, cons_temp, cmap = :inferno)
plt.colorbar(mappable = p, label = "Conservative Temperature (ºC)")
gca().invert_yaxis()
plt.title("Conservative Temperature with Mixed Layer Line")
plt.xlabel("Time")
plt.ylabel("Depth(m)")

#Handle time ticks
t_ticks = range(extrema(times)...; step=(maximum(times) - minimum(times))/8)
t_tickvals = Dates.value.(t_ticks)
t_ticklabels = Dates.format.(t_ticks, dateformat"yyyy-mm-dd")
plt.xticks(ticks=t_tickvals, labels = t_ticklabels, rotation=20)

plt.plot(Dates.value.(times)', mixed_depths; label="Mixed Layer 0.03kg/m^3 threshold")

plt.tight_layout()
savefig(joinpath(visdir, "mixed_layer_depth_on_temp.png"))
plt.close()

p = contourf(num_times, depths, abs_sal, cmap = "cividis")
plt.colorbar(mappable = p, label = "Absolute Salinity (g/kg)")
gca().invert_yaxis()
plt.title("Absolute Salinity with Mixed Layer Line")
plt.xlabel("Time")
plt.ylabel("Depth(m)")

#Handle time ticks
t_ticks = range(extrema(times)...; step=(maximum(times) - minimum(times))/8)
t_tickvals = Dates.value.(t_ticks)
t_ticklabels = Dates.format.(t_ticks, dateformat"yyyy-mm-dd")
plt.xticks(ticks=t_tickvals, labels = t_ticklabels, rotation=20)

plt.plot(Dates.value.(times)', mixed_depths; label="Mixed Layer 0.03kg/m^3 threshold")

plt.tight_layout()
savefig(joinpath(visdir, "mixed_layer_depth_on_sal.png"))
plt.close()