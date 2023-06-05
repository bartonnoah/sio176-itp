using LinearAlgebra, MAT, Plots, StatsBase

cd(@__DIR__)
datafile = "../data/interpolated_vars.mat"

#Load in the data
data = matread(datafile)
cons_temp = data["cons_temp"]
abs_sal = data["abs_sal"]
presvals = data["pres"][:,1] #pressure is a regular grid

#Do a PCA using conservative temp and abs salinity as the two variables of interest
#Flatten and combine the datasets, and then center and whiten
pca_data = hcat(cons_temp', abs_sal')
dataset_indicators = hcat((fill(i, size(cons_temp')) for i in 1:2)...)
pca_data .-= mean(pca_data, dims = 1)

#whiten, accounting for the fact that small temp variations shouldn't count the same as large ones
scale_dict = Dict()
for i in 1:2
    subset = @view pca_data[dataset_indicators .== i]
    subset ./= std(subset)
    scale_dict[i] = std(subset)
end

#Calculate the svd
svd_obj = svd(pca_data; full=false)
U, S, V = svd_obj.U, svd_obj.S, svd_obj.V

perc_var = round.(Int, 100 .* S ./ sum(S))

#now get the principal components
principial_components = pca_data * V

for pcnum in 1:5
    pc = V[:, pcnum]
    tempsize = size(cons_temp, 1)
    pctemp = pc[1:tempsize] * scale_dict[1]
    pcsal = pc[tempsize+1:end] * scale_dict[2]
    display(pctemp)
    display(pcsal)

    #Now plot the data
    p1 = plot(pctemp, presvals; yflip=true, label="")
    plot!(p1, xlabel = "Conservative Temp Deviation(ºC) ", ylabel = "Pressure (dbar)")
    p2 = plot(pcsal, presvals; yflip=true, label="") 
    plot!(p2, xlabel = "Absolute Salinity Deviation (g/kg) ", ylabel = "Pressure (dbar)")

    #Now put them into a layout
    l = @layout [a b]
    p = plot(p1, p2, layout = l, plot_title = "Principal Component $pcnum ($(perc_var[pcnum])% of variance)", plot_title_fontsize = 12)
    savefig("PC$(pcnum).png")
end


