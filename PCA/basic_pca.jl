using LinearAlgebra, MAT, Plots, StatsBase, Dates

gr()
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
pca_datasets = (cons_temp', abs_sal')
pca_data = hcat(pca_datasets...)
dataset_indicators = hcat((fill(i, size(dataset)) for (i, dataset) in enumerate(pca_datasets))...)
subsetnums = 1:length(pca_datasets)
pca_data .-= mean(pca_data, dims = 1)

#whiten, accounting for the fact that small temp variations shouldn't count the same as large ones
scalefactor = []
for i in subsetnums
    subset = @view pca_data[dataset_indicators .== i]
    dev = std(subset)
    subset ./= dev
    push!(scalefactor, dev)
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
    pctemp = pc[1:tempsize] * scalefactor[1]
    pcsal = pc[tempsize+1:2*tempsize] * scalefactor[2]

    #Now plot the data
    p1 = plot(pctemp, presvals; yflip=true, label="", rotation=20)
    plot!(p1, xlabel = "Conservative Temp Deviation(ºC) ", ylabel = "Pressure (dbar)")
    p2 = plot(pcsal, presvals; yflip=true, label="", rotation=20) 
    plot!(p2, xlabel = "Absolute Salinity Deviation (g/kg) ")

    #Now put them into a layout
    l = @layout [a b]
    p = plot(p1, p2, layout = l, plot_title = "Principal component $pcnum ($(perc_var[pcnum])% of variance)", plot_title_fontsize = 12, figsize = (10,10))
    savefig(p, "PC$(pcnum).png")
end

#Now plot a time series of the PCA components
times = @. Millisecond(round(Int, data["time"] * milliseconds_in_day)) + DateTime(0,1,1)

for i in 1:3
    plot(;title="PC $i over time", xlabel = "Time", ylabel = "Value")
    plot!(times[:], principal_components[:, i], label="")
    savefig("pca$(i)_over_time.png")
end




