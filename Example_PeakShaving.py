# -*- coding: utf-8 -*-
"""
@author: Joel Alpízar Castillo
Version: 1.0
Update: 07-09-22
"""

from PeakShaving import GaussPV, BESS, PeakShaving, met, load, NREL_PV_model, PV_system, year2months, yearlycosts
import csvreader


from numpy import array, arange
import matplotlib.pyplot as plt
import time

###############################################################################    

start = time.time()    # The timer is initializad.
## Data from a case scenario in Costa Rica
SpecificEnergy = array([4.97, 5.3, 5.24, 4.39, 3.81, 3.67, 3.67, 3.81, 3.87, 3.67, 3.74, 4.31])
SunriseTime = array([6.02, 6.00, 5.46, 5.27, 5.16, 5.17, 5.24, 5.29, 5.28, 5.28, 5.35, 5.50])
SunsetTime = array([17.37, 17.49, 17.52, 17.52, 17.57, 18.05, 18.08, 18.00, 17.42, 17.23, 17.14, 17.21])
NoonTime = array([12.45, 12.4, 12.4, 12.3, 12.3, 12.3, 12.3, 12.25, 12.05, 11.5, 12, 12.2])
Deviations = array([2.72, 2.95, 2.85, 2.95, 2.6, 2.6, 2.8, 2.8, 2.65, 2.55, 2.6, 2.6])
Alpha = array([-3, 3])
DCPeakPower = 17*0.275 # Optimizers

# Load measurements
#CSVData = csvreader.read_data(csv='L_EE_5kWh.csv', address='')
#CSVData = csvreader.read_data(csv='L_EE_10kWh.csv', address='')
#CSVData = csvreader.read_data(csv='L_EE_20kWh.csv', address='')
#CSVData = csvreader.read_data(csv='L_EE_30kWh.csv', address='')
CSVData = csvreader.read_data(csv='L_EE_40kWh.csv', address='')
#CSVData = csvreader.read_data(csv='L_EE_50kWh.csv', address='')

# PV measurement
CSVData_PV = csvreader.read_data(csv='UCR_PV_measurements.csv', address='')


# BESS data

###############################################################################
## Main code excecution
###############################################################################


#met_matrix='EIE_1h.mat'
met_matrix='EIE_2020_1h.mat'
#load_matrix='L_CR.mat'

met_data = met()
met_data.load_data(met_matrix)
#load_data = load()
#load_data.load_data(load_matrix)

      
 ## UCR ##     
 # Given by the design.        
Opt_theta = 10.5 # tilt of the module
Opt_Am = 200 # azimuth of the module
alpha = 0.2 # This is normally the Albedo coefficient used, it can be changed based on PV location  

 # Given by PV module datasheet
eff_datasheet = 0.1619 #given by the datasheet
C_eff = -0.0035 # thermal coeficient 
thao_alpha = 0.9
Tm_NOCT = 45 #Temperature of the PV module at NOCT, data sheet  
Ta_NOCT = 20
Voc_STC = 38.6
Isc_STC = 9.03
Vmpp_STC = 31.4
Impp_STC = 8.44
Area = 1.65*0.992

# Given by the inverter/charger datasheet
eff_mppt = 0.96
eff_inverter = 0.94

# Given by the battery datasheet
V_nom = 3.6 #nominal voltage of the electrochemical cell Panasonic NCR18650GA  
n_serie = 125 # number of cells in series
cell_cap = 3.45 # Capacity of the cell
n_parallel = 5 # number of cells in parallel
n_batteries = 5 # number of batteries
Capacity = 10.78125
#Capacity = cell_cap*n_serie*n_parallel*n_batteries*V_nom/1000 # Capacity of the battery in kWh
eff_battery = 0.97; # battery efficiency
SoC_lowerlimit = 0.2 # lower limit to define the range of operation of the battery
SoC_upperlimit = 0.8 # upper limit to define the range of operation of the battery



PV_system = PV_system()
PV_system.load_module_orientation(Opt_theta, Opt_Am, alpha)
PV_system.load_module_datasheet(eff_datasheet, C_eff, thao_alpha, Tm_NOCT, Ta_NOCT, Voc_STC, Isc_STC, Vmpp_STC, Impp_STC, Area)
PV_system.load_inverter_datasheet(eff_mppt, eff_inverter)
PV_system.load_battery_datasheet(V_nom, n_serie, cell_cap, n_parallel, eff_battery, SoC_lowerlimit, SoC_upperlimit)


# Battery at UCR
Battery = BESS(Capacity, MaxPower = 0.5, Char_eff = eff_battery, Dischar_eff = eff_battery, Self_discharge = 0, MinSoC = SoC_lowerlimit, MaxSoC = SoC_upperlimit)

# Time variables

step = 1; # this is because the simulation is running with 1h timestep
Time_zone = 'US/Central'
# ## UCR ##
# # First the longitude and latitude are needed as well as the date and time
L=9.936926 #this is the latitude of the observer
l=-84.0439 #this is the longitude of the observer
N_m=17 #number of PV modules of the PV system 

# Simulation starting date 
year_s = 2020
month_s = 1
day_s = 1
hour_s = 1
minute_s = 0
second_s = 0

# Simulation finishing date 
year_e = 2020
month_e = 1
day_e = 1
hour_e = 0
minute_e = 0
second_e = 0


PV_model = NREL_PV_model()
PV_model.sim_init(met_data, PV_system, step, year_s, month_s, day_s, hour_s, minute_s, second_s, year_e, month_e, day_e, hour_e, minute_e, second_e, L, l, N_m, Time_zone)

#end = time.time()  # The simulation time is saved.
#totalelapsed = end - start  # The total time is calculated.

#PV_out = PV_model.P_pv_syst
#fig=plt.figure(5)
#plt.plot(PV_out)
#
#plt.ylabel('PV Power, P [W]')
#plt.xlabel('Time, t [h]')
#plt.grid(True)
#plt.show()
###############################################################################


SymData = GaussPV(SpecificEnergy, SunriseTime, SunsetTime, NoonTime, Deviations, Alpha, TimeStep = 0.25)
SymData.EstimateDCPowerOutput(DCPeakPower)

PermPower = 1


CSVData.data2array()
CSVData_PV.data2array()
## EE
LoadData0 = CSVData.ar[:,0]/1000
LoadData1 = [0]*len(LoadData0)*4


for i in range(len(LoadData0)):
    for j in range(4):
        LoadData1[4*i+j] = LoadData0[i]
        

## UCR
PVMeasurements = CSVData_PV.ar[:,0]

Letter_size = 16

plt.rc('font', size=Letter_size)          # controls default text sizes
plt.rc('axes', titlesize=Letter_size)     # fontsize of the axes title
plt.rc('axes', labelsize=Letter_size)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=Letter_size)    # fontsize of the tick labels
plt.rc('ytick', labelsize=Letter_size)    # fontsize of the tick labels
plt.rc('legend', fontsize=Letter_size)    # legend fontsize
plt.rc('figure', titlesize=Letter_size)  # fontsize of the figure title


SimulationResults1 = PeakShaving(PVMeasurements[0:4*24*365], Battery, LoadData1[0:4*24*365], MaxPower = PermPower, EndTime = 24*365, TimeStep = 0.25)
SimulationResults1.Initialize(InitialSoC = 0.5)
SimulationResults1.Simulate()

SoC1 = SimulationResults1.BESS.SoC
#
fig1 = plt.figure(1)
plt.plot(SimulationResults1.Load[0:4*24*3], 'r', label='Load')
plt.plot(SimulationResults1.PV[0:4*24*3], 'b', label='PV')
plt.plot(SimulationResults1.BESS.Power[0:4*24*3], 'k', label='BESS')
plt.plot(SimulationResults1.Grid[0:4*24*3], 'g', label='Grid')
plt.xlim([0, 4*24*3])
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
plt.xlabel('Days')
#plt.yticks([-1, 0, 1, 1.5, 2, 3, 4])
plt.yticks([-1, 0, 1, 2, 3, 4])
plt.ylim([-1.5, 4])
plt.ylabel('Power [kW]')
#plt.title('Power behaviour with real PV measurements')
plt.legend(loc='upper left', ncol=4)
plt.grid()
plt.show()

SimulationResults2 = PeakShaving(SymData.Power[0:4*24*365], Battery, LoadData1[0:4*24*365], MaxPower = PermPower, EndTime = 24*365, TimeStep = 0.25)
SimulationResults2.Initialize(InitialSoC = 0.5)
SimulationResults2.Simulate()

SoC2 = SimulationResults2.BESS.SoC

fig2 = plt.figure(2)
plt.plot(SimulationResults2.Load[0:4*24*3], 'r', label='Load')
plt.plot(SimulationResults2.PV[0:4*24*3], 'b', label='PV')
plt.plot(SimulationResults2.BESS.Power[0:4*24*3], 'k', label='BESS')
plt.plot(SimulationResults2.Grid[0:4*24*3], 'g', label='Grid')
plt.xlim([0, 4*24*3])
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
plt.xlabel('Days')
#plt.yticks([-1, 0, 1, 1.5, 2, 3, 4])
plt.yticks([-1, 0, 1, 2, 3, 4])
plt.ylim([-1.5, 4])
plt.ylabel('Power [kW]')
#plt.title('Power behaviour with the Gaussian model')
plt.legend(loc='upper left', ncol=4)
plt.grid()
plt.show()

SimulationResults3 = PeakShaving(array(PV_model.P_pv_syst[0:4*24*365]), Battery, LoadData1[0:4*24*365], MaxPower = PermPower, EndTime = 24*365, TimeStep = 0.25)
SimulationResults3.Initialize(InitialSoC = 0.5)
SimulationResults3.Simulate()

SoC3 = SimulationResults3.BESS.SoC

fig3 = plt.figure(3)
plt.plot(SimulationResults3.Load[0:4*24*3], 'r', label='Load')
plt.plot(SimulationResults3.PV[0:4*24*3], 'b', label='PV')
plt.plot(SimulationResults3.BESS.Power[0:4*24*3], 'k', label='BESS')
plt.plot(SimulationResults3.Grid[0:4*24*3], 'g', label='Grid')
plt.xlim([0, 4*24*3])
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
plt.xlabel('Days')
#plt.yticks([-1, 0, 1, 1.5, 2, 3, 4])
plt.yticks([-1, 0, 1, 2, 3, 4])
plt.ylim([-1.5, 4])
plt.ylabel('Power [kW]')
#plt.title('Power behaviour with the NREL model')
plt.legend(loc='upper left', ncol=4)
plt.grid()
plt.show()


fig4 = plt.figure(4)
plt.plot(SoC1[0:4*24*3], 'r', label='PV measurements')
plt.plot(SoC2[0:4*24*3], 'b', label='Gaussian model')
plt.plot(SoC3[0:4*24*3], 'k', label='NREL model')
plt.xlim([0, 4*24*3])
plt.xlabel('Days')
plt.ylabel('SoC [%]')
#plt.title('SoC of the BESS')
plt.legend(loc='upper left', ncol=4)
plt.grid(axis='y')
plt.show()


e1 = SoC1 - SoC2
e2 = SoC1 - SoC3

fig5 = plt.figure(5)
plt.plot(e1, 'r', label='Error with the Gaussian model') # , lhttps://www.msn.com/en-us/feed?lacc=11abel=''
plt.plot(e2, 'b', label='Error with the NREL model') # , label=''
plt.xlim([0, 365*24*4])
plt.xlabel('Time steps')
plt.ylabel('SoC [%]')
#plt.title('Difference in SoC')
plt.legend(loc='upper right')
plt.grid()
plt.show()


fig6 = plt.figure(6)
plt.plot(SimulationResults1.PV[4*24*177:4*24*184], 'r', label='PV measurements')
plt.plot(SimulationResults2.PV[4*24*177:4*24*184], 'b', label='Gaussian model')
plt.plot(SimulationResults3.PV[4*24*177:4*24*184], 'k', label='NREL model')
plt.xlim([0, 4*24*7])
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
plt.xlabel('Days')
#plt.xticks(arange(7), ['06/26', '06/27', '06/28', '06/29', '06/30', '07/01', '07/02'])
plt.ylabel('PV Power [kW]')
#plt.title('monthly Costs')
plt.legend(loc='upper left', ncol=3)
plt.grid(axis='y')
plt.show()

###############################################################################

dt = 0.25

CRC2USD = 570

cost_hours = [6, 10.5, 12.5, 17, 20, 24]
cost_night = 26.48 / CRC2USD
cost_valley = 63.28 / CRC2USD
cost_peak = 154.35 / CRC2USD
costs_energy = [cost_night, cost_valley, cost_peak, cost_valley, cost_peak, cost_night]
#costs_power = [0 for i in costs_energy]
costs_power = 0

year_cost_list = yearlycosts(SimulationResults1.Grid, cost_hours, costs_energy)


monthly_power_Measurements = year2months(SimulationResults1.Grid, dt)
monthly_power_cost_Measurements = [max(month)*costs_power for month in monthly_power_Measurements]
year_energy_cost_Measurements = [dt*SimulationResults1.Grid[i]*year_cost_list[i] for i in range(len(year_cost_list))]
month_energy_cost_Measurements = year2months(year_energy_cost_Measurements, dt)
monthly_energy_cost_Measurements = [sum(month) for month in month_energy_cost_Measurements]
monthly_energy_Measurements = [dt*sum(month) for month in year2months(SimulationResults1.PV, dt)]

monthly_power_Gauss = year2months(SimulationResults2.Grid, dt)
monthly_power_cost_Gauss = [max(month)*costs_power for month in monthly_power_Gauss]
year_energy_cost_Gauss = [dt*SimulationResults2.Grid[i]*year_cost_list[i] for i in range(len(year_cost_list))]
month_energy_cost_Gauss = year2months(year_energy_cost_Gauss, dt)
monthly_energy_cost_Gauss = [sum(month) for month in month_energy_cost_Gauss]
monthly_energy_Gauss = [dt*sum(month) for month in year2months(SimulationResults2.PV, dt)]

monthly_power_NREL = year2months(SimulationResults3.Grid, dt)
monthly_power_cost_NREL = [max(month)*costs_power for month in monthly_power_NREL]
year_energy_cost_NREL = [dt*SimulationResults3.Grid[i]*year_cost_list[i] for i in range(len(year_cost_list))]
month_energy_cost_NREL = year2months(year_energy_cost_NREL, dt)
monthly_energy_cost_NREL = [sum(month) for month in month_energy_cost_NREL]
monthly_energy_NREL = [dt*sum(month) for month in year2months(SimulationResults3.PV, dt)]

fig7 = plt.figure(7)
plt.plot(monthly_energy_cost_Measurements, 'r', label='PV measurements')
plt.plot(monthly_energy_cost_Gauss, 'b', label='Gaussian model')
plt.plot(monthly_energy_cost_NREL, 'k', label='NREL model')
plt.xlim([0, 11])
plt.xlabel('Months')
plt.xticks(arange(12), ['Jan', 'Feb', 'March', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
plt.ylabel('Costs [CRC]')
plt.title('Monthly Costs')
plt.legend(loc='lower center')
plt.grid()
plt.show()

#monthly_energy = [monthly_energy_Measurements, monthly_energy_Gauss, monthly_energy_NREL]

#X = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
#X_axis = arange(len(X))
#fig8 = plt.figure(8)
#plt.bar(X_axis - 0.2, monthly_energy_Measurements, 0.2, color = 'r', label = 'PV measurements')
#plt.bar(X_axis + 0.0, monthly_energy_Gauss, 0.2, color = 'b', label = 'Gaussian model')
#plt.bar(X_axis + 0.2, monthly_energy_NREL, 0.2, color = 'k', label = 'MDB model')
#plt.xlabel('Months')
#plt.xticks(X_axis, X)
##plt.xticks(arange(12), ['Jan', 'Feb', 'March', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
#plt.ylabel('Energy [kWh]')
#plt.title('Monthly PV generation')
#plt.legend(loc='upper center', ncol=3)
#plt.grid()
#plt.show()

monthly_power_Load = year2months(SimulationResults3.Load, dt)
monthly_power_cost_Load = [max(month)*costs_power for month in monthly_power_Load]
year_energy_cost_Load = [dt*SimulationResults3.Load[i]*year_cost_list[i] for i in range(len(year_cost_list))]
month_energy_cost_Load = year2months(year_energy_cost_Load, dt)
monthly_energy_cost_Load = [sum(month) for month in month_energy_cost_Load]





end = time.time()  # The simulation time is saved.
totalelapsed = end - start  # The total time is calculated.