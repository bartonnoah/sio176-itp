using MAT, GibbsSeaWater, PyCall, StaticArrays, StatsBase
@pyimport scipy.interpolate as spitp

cd(@__DIR__)

data = matread("itp121_data.mat")

times = data["time"]
nice_times = range(extrema(times)...; length = length(times))[:]'

pres = data["pres"]

#Expand times into a matrix
expanded_tvals = repeat(times, size(pres, 1))

#Interpolate pressure (it's supposed to be a regular grid but some values are marked NaN for some reason)
pres[:, :] .= 6.0:2.:760.
nice_pres = 6.0:2.:760.

#Make the interpolation grids
interp_grid = hcat(vec(pres), vec(expanded_tvals))
original_shape = size(pres)

#Interpolate temp
temp = data["temperature"]
temp_na = vec(isnan.(temp))
temp_not_na = (!).(temp_na)

temp_itp = spitp.griddata(interp_grid[temp_not_na, :], temp[temp_not_na], interp_grid[temp_na, :])
temp[temp_na] .= temp_itp

#Interpolate psu
sal = data["salinity"]
sal[sal.>=2000] .= NaN
sal_na = vec(isnan.(sal))
sal_not_na = (!).(sal_na)

sal_itp = spitp.griddata(interp_grid[sal_not_na, :], sal[sal_not_na], interp_grid[sal_na, :])
sal[sal_na] .= sal_itp

#Now trim the arrays to remove the NaNs that persist on the edges of the dataset due to being outside the convex hull
idx = (Colon(), 2:size(times, 2)-1)
time_trimmed = times[idx...]
pres_trimmed = pres[idx...]
temp_trimmed = temp[idx...]
sal_trimmed = sal[idx...]
lat_trimmed = data["latitude"][idx...]
lon_trimmed = data["longitude"][idx...]

#Now generate the derived derived variables
abs_sal = gsw_sa_from_sp.(sal_trimmed, pres_trimmed, lon_trimmed, lat_trimmed)
cons_temp = gsw_ct_from_t.(abs_sal, temp_trimmed, pres_trimmed)
depth = -1 .* gsw_z_from_p.(pres_trimmed, lat_trimmed, 0, 0)
pot_dens_anom = gsw_sigma0.(abs_sal, cons_temp)

#Now calculate mixed layer depth
diff_array = copy(pot_dens_anom)
mixed_pres = permutedims(zeros(size(diff_array, 2))')

for (i,col) in enumerate(eachcol(diff_array))
    col .-= col[5]
    mixed_idx = findfirst(x -> (x>0.03), col)
    if (mixed_idx) === nothing || (findfirst(!isnan, col) > 10)
        mixed_depths[i] = NaN
        col .= NaN
        continue
    end
    mixed_pres[i] = pres_trimmed[mixed_idx, i]
end

datadict = Dict("time"=>time_trimmed, "pres"=>pres_trimmed, "temp"=>temp_trimmed, "sal"=>sal_trimmed, 
                "lat"=>lat_trimmed, "lon"=>lon_trimmed, "abs_sal"=>abs_sal, "cons_temp"=>cons_temp, "depth"=>depth,
                "pot_dens_anom"=>pot_dens_anom, "Description"=>"Interpolated values for all major variables, with first and last profile dropped due to inability to interpolate them properly", "mixed_layer_pressure"=>mixed_pres )

matwrite("interpolated_vars.mat", datadict)


