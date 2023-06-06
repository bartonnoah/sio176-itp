using LinearAlgebra, MAT, Plots, StatsBase, Dates

cd(@__DIR__)
datafile = "../data/interpolated_vars.mat"

#Load in the data
data = matread(datafile)
cons_temp = data["cons_temp_denoised"]
abs_sal = data["abs_sal_denoised"]
presvals = data["pres"][:,1] #pressure is a regular grid

milliseconds_in_day = Dates.value(Millisecond(Day(1)))
times = @. Millisecond(round(Int, data["time"] * milliseconds_in_day)) + DateTime(0,1,1)

#Do a PCA using conservative temp and abs salinity as the two variables of interest
#Flatten and combine the datasets, and then center and whiten
pca_datasets = (cons_temp', abs_sal', data["mixed_layer_pressure"])
pca_data = hcat(pca_datasets...)
dataset_indicators = hcat((fill(i, size(dataset)) for (i, dataset) in enumerate(pca_datasets))...)
subsetnums = 1:length(pca_datasets)
pca_data .-= mean(pca_data, dims = 1)

#whiten, accounting for the fact that small temp variations shouldn't count the same as large ones
for i in subsetnums
    subset = @view pca_data[dataset_indicators .== i]
    subset ./= mad(subset; center=0.)
end

#Calculate the svd
svd_obj = svd(pca_data; full=false)
U, S, V = svd_obj.U, svd_obj.S, svd_obj.V

perc_var = round.(Int, 100 .* S ./ sum(S))

#now get the principal components
principal_components = pca_data * V

for pcnum in 1:5
    pc = V[:, pcnum]
    tempsize = size(cons_temp, 1)
    scalefactor = 25
    pctemp = pc[1:tempsize] * scalefactor
    pcsal = pc[tempsize+1:2*tempsize] * scalefactor
    pcmixed = pc[end] * scalefactor

    #Now plot the data
    p1 = plot(pctemp, presvals; yflip=true, label="")
    plot!(p1, xlabel = "Conservative Temp Deviation(ºC) ", ylabel = "Pressure (dbar)")
    p2 = plot(pcsal, presvals; yflip=true, label="") 
    annotate!(p2, ((0.4, 0.9), Plots.text("Mixed Layer: $(round(pcmixed, digits=4))dbar", 10)))
    plot!(p2, xlabel = "Absolute Salinity Deviation (g/kg) ")

    #Now put them into a layout
    l = @layout [a b]
    p = plot(p1, p2, layout = l, plot_title = "Principal component $pcnum ($(perc_var[pcnum])% of variance)", plot_title_fontsize = 12, figsize = (10,10))
    savefig(p, "PC$(pcnum).png")
end

#Now plot a time series of the PCA components
times = @. Millisecond(round(Int, data["time"] * milliseconds_in_day)) + DateTime(0,1,1)

myp = plot(;title="PCA Components over time", xlabel = "Time", ylabel = "Value")
for i in 1:5
    plot!(times[:], principal_components[:, i]; label = "$i")
end

savefig(myp, "pca_over_time.png")


