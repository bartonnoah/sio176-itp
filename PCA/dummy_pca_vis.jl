using LinearAlgebra, MAT, Plots, StatsBase, Dates

pyplot()
cd(@__DIR__)
datafile = "../data/interpolated_vars.mat"

#Load in the data
data = matread(datafile)
cons_temp = data["cons_temp_denoised"]
abs_sal = data["abs_sal_denoised"]

milliseconds_in_day = Dates.value(Millisecond(Day(1)))
times = @. Millisecond(round(Int, data["time"] * milliseconds_in_day)) + DateTime(0,1,1)

#Do a PCA using conservative temp and abs salinity as the two variables of interest
#Flatten and combine the datasets, and then center and whiten
decibars = 60
idx = decibars ÷ 2 - 2
saldat = abs_sal[idx,:]
tempdat = cons_temp[idx,:]
pca_datasets = (tempdat, saldat)
pca_data = hcat(pca_datasets...)
dataset_indicators = hcat((fill(i, size(dataset)) for (i, dataset) in enumerate(pca_datasets))...)
subsetnums = 1:length(pca_datasets)
data_means = mean(pca_data, dims = 1)
pca_data .-= data_means

#Calculate the svd
svd_obj = svd(pca_data; full=false)
U, S, V = svd_obj.U, svd_obj.S, svd_obj.V

#perc_var = round.(Int, 100 .* S ./ sum(S))

#now get the principal components
principal_components = pca_data * V

deviations_from_mean(x) = (x.-mean(x))./std(x)

#Now plot original data
p = scatter(tempdat, saldat; xlabel = "Conservative Temperature (ºC)", ylabel = "Absolute Salinity (g/kg)", 
title = "Sample PCA: Temp and Sal at $(decibars)db", label = "", aspect_ratio = :equal)
quiver!(p, [data_means[1] ], [data_means[2] ], quiver = (V[1, [1]], V[2, [1]]); lw = 4, label = "PC 1")
quiver!(p, [data_means[1] ], [data_means[2] ], quiver = (V[1, [2]], V[2, [2]]); lw = 4, label = "PC 2")

savefig(p, "example_pca.png")



