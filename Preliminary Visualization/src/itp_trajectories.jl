using MAT, PyCall, GibbsSeaWater, Dates, ColorSchemes, StatsBase
@pyimport cartopy.crs as ccrs
@pyimport matplotlib.pyplot as plt
@pyimport cartopy.mpl.ticker as cmplt
@pyimport cartopy.mpl.gridliner as cmplg
@pyimport matplotlib.ticker as mticker

LongitudeFormatter, LatitudeFormatter = cmplt.LongitudeFormatter, cmplt.LatitudeFormatter
LONGITUDE_FORMATTER, LATITUDE_FORMATTER = cmplg.LONGITUDE_FORMATTER, cmplg.LATITUDE_FORMATTER

cd(@__DIR__)

data = matread("../../data/itp121_data.mat")

times = @. Millisecond(round(Int, data["time"] * Dates.value(Millisecond(Day(1))))) + DateTime(0,1,1)

lats = data["latitude"]

lons = data["longitude"]

proj = ccrs.NorthPolarStereo()

ax = plt.subplot(111, projection = proj)
ax.coastlines()

mappable = ax.scatter(lons, lats; c = Dates.value.(times), transform = ccrs.PlateCarree(), s=1)

tick_dates = range(extrema(times)...; step = Month(3))
tick_vals = Dates.value.(tick_dates)
tick_labels = Dates.format.(tick_dates, dateformat"yyyy-mm-dd")

cbar = plt.colorbar(mappable; ticks = tick_vals)
cbar.ax.set_yticklabels(tick_labels)


plt.scatter([lons[begin]],[lats[begin]]; label="Start", color="pink", transform = ccrs.PlateCarree(),s=2)
plt.scatter([lons[end]],[lats[end]]; label="End", color="red", transform = ccrs.PlateCarree(),s=2)
plt.gca().gridlines(draw_labels=true, dms=true, x_inline=false, y_inline=false)
plt.legend()
loncirc = -180:2:180
lat = fill(60, length(loncirc))
ax.scatter(loncirc, lat; transform=ccrs.PlateCarree(), alpha=0)

plt.title("Drifter Positions as a function of time")
plt.tight_layout()

plt.savefig("../vis/drifter_path.png")

plt.close()
