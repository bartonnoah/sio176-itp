using MAT, GibbsSeaWater, PyCall, StaticArrays, StatsBase, ImageFiltering
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

temp_itp = spitp.griddata(interp_grid[temp_not_na, :], temp[temp_not_na], interp_grid[temp_na, :]; method="cubic",rescale = true)
temp[temp_na] .= temp_itp

#Interpolate psu
sal = data["salinity"]
sal[sal.>=2000] .= NaN
sal_na = vec(isnan.(sal))
sal_not_na = (!).(sal_na)

sal_itp = spitp.griddata(interp_grid[sal_not_na, :], sal[sal_not_na], interp_grid[sal_na, :]; method = "cubic", rescale = true)
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
    col .-= first(col)
    mixed_idx = findfirst(x -> (x>0.03), col)
    mixed_pres[i] = pres_trimmed[mixed_idx, i]
end

#Now finally add a median filter version to reduce noise
patch_size = (9,21)
denoised_cons_temp = mapwindow(median, cons_temp, patch_size)
denoised_abs_sal =mapwindow(median, abs_sal, patch_size)
denoised_pot_dens =mapwindow(median, pot_dens_anom, patch_size)

#Now calculate mixed layer depth
smooth_diffs = copy(denoised_pot_dens)
smooth_mixed_pres = permutedims(zeros(size(diff_array, 2))')

for (i,col) in enumerate(eachcol(smooth_diffs))
    col .-= first(col)
    mixed_idx = findfirst(x -> (x>0.03), col)
    smooth_mixed_pres[i] = pres_trimmed[mixed_idx, i]
end


datadict = Dict("time"=>time_trimmed, "pres"=>pres_trimmed, "temp"=>temp_trimmed, "sal"=>sal_trimmed, 
                "lat"=>lat_trimmed, "lon"=>lon_trimmed, "abs_sal"=>abs_sal, "cons_temp"=>cons_temp, "depth"=>depth,
                "pot_dens_anom"=>pot_dens_anom, "Description"=>"Interpolated values for all major variables, with first and last profile dropped due to inability to interpolate them properly", "mixed_layer_pressure"=>mixed_pres,
                "cons_temp_denoised"=>denoised_cons_temp, "abs_sal_denoised"=>denoised_abs_sal, "pot_dens_anom_denoised"=>denoised_pot_dens,
                "smooth_mixed_layer_pressure"=>smooth_mixed_pres)

matwrite("interpolated_vars.mat", datadict)


