# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 14:45:45 2017

@author: zagoven
"""
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
import matplotlib.pyplot as plt
import numpy as np
#from datetime 
import datetime 

#choose the year,month and day in our netcdf
year=2014
month=11
day=24
# input netcdf:
filein='ROMS_Laptev_Sea_NETCDF3_CLASSIC_each_day.nc'
# output file for using in methan.f90 model:
fileout='Laptev_sea_{0}_{1}_{2}'.format(day,month, year)+'.dat'
# output temperature plot:
plot_T='Temp_Laptev_sea_{0}_{1}_{2}'.format(day,month, year)+'.png' 
# output salinity plot:
plot_S='Sal_Laptev_sea_{0}_{1}_{2}'.format(day,month, year)+'.png'
# create the human formate for our time array
dt=datetime.datetime(year,month,day,12,0)
#reading the netcdf
f = Dataset('ROMS_Laptev_Sea_NETCDF3_CLASSIC_each_day.nc')
depth =  f.variables['depth'][:] # means all value
sal = f.variables['sal'][:]
temp = f.variables['temp'][:]
ocean_time =  f.variables['time'][:]
times = num2date(ocean_time,
                 units= 'seconds since 1948-01-01 00:00:00' ,
                                   calendar= 'standard')

#np.where(times==dt) #ищет в массивах, где наше время совпадает с dt
#sort the depth for methan.f90 
#(not sure that nessesary - need to check)
True_depth = -np.sort(depth)
ind=np.where(times==dt)[0][0] #index in time, [0] for the acsess to the 3 element

sal[ind] #salinity from one particular time moment
temp[ind] # the same for temperature
#create 3 columns from 3 arrays
out_array=np.column_stack((True_depth,sal[ind],temp[ind])) 

np.savetxt(fileout,out_array, delimiter=' ') #save the main data file

plt.plot(temp[ind],depth) #creation temperature plot
plt.ylim(150,0)
plt.title(str(times[ind])+' Temperature')
plt.xlabel('Temperature,°С')
plt.ylabel('Depth,m')
 #it is nessesary to put .savefig before .show
plt.savefig(plot_T, dpi=400, facecolor='w', edgecolor='r',
        orientation='portrait', papertype=None, forma=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)
plt.show()

#creation salinity plot
plt.plot(sal[ind],depth)
plt.ylim(150,0)
plt.title(str(times[ind])+' Salinity')
plt.xlabel('Salinity,psu')
plt.ylabel('Depth,m')
plt.savefig(plot_S, dpi=400, facecolor='w', edgecolor='r',# doesn't draw the line in the export plot.png
        orientation='portrait', papertype=None, forma=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)
plt.show()


