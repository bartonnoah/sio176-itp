using MAT, Plots, GibbsSeaWater, Dates, ColorSchemes

pyplot()
cd(@__DIR__)

datafile = "../../data/itp121_data.mat"
derived_file = "../../data/derived_variables.mat"
visdir = "../vis"

raw_vars = matread(datafile)
derived_vars = matread(derived_file)

derived_vars["abs_sal"][derived_vars["abs_sal"] .> 2000] .= NaN

times = @. Millisecond(round(Int, raw_vars["time"] * Dates.value(Millisecond(Day(1))))) + DateTime(0,1,1)

comp_dates = (times[320], times[700])

comp_names = Dates.format.(comp_dates, dateformat"yyyy-mm-dd")

comp_idxs = findfirst.([==(d) for d in comp_dates], Ref(times[:]))

for (var, varname, unit) in zip(["cons_temp", "abs_sal"], ["Conservative Temperature", "Absolute Salinity"], ["ºC", "g/kg"])
    myp = plot(;title = "$varname on $(comp_names[1]) and $(comp_names[2])", ylabel = "Depth (m)", xlabel = "$varname ($unit)")
    for (idx, dayname) in zip(comp_idxs, comp_names)
        plot!(myp, derived_vars[var][:, idx], derived_vars["depth"][:, idx]; label = dayname, yflip = true, dpi=1200)
    end
    savefig(myp, joinpath(visdir, "$varname comparison.png"))
end