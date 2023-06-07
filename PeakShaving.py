# -*- coding: utf-8 -*-
"""
@author: Joel Alpízar Castillo
Version: 1.0
Update: 23-08-22
"""


#import csvreader

from numpy import array
import scipy.io
import sympy as sp


class GaussPV:
    def __init__(self, SpecificEnergy, SunriseTime, SunsetTime, NoonTime, Deviations, Alpha, EndTime=24, TimeStep = 0.25, StartTime = 0):
        from numpy import arange
        # TimeStep is the fraction of hour used to create the TimeStamp
        self.TimeStamp = arange(StartTime, EndTime, TimeStep)              
        self.SunriseTime = SunriseTime
        self.SunsetTime = SunsetTime
        self.NoonTime = NoonTime
        self.Deviations = Deviations
        self.Alpha = Alpha
        self.SpecificEnergy = SpecificEnergy
        
    def Transform2DecimalHour(self, TimeList):
        return array([element//1 + (element%1)/0.6 for element in TimeList])

    def GaussPVEnvelopes(self):
        from numpy import sqrt, exp, logical_and, pi # , erf
        from scipy.special import erf
                
        return array([[((self.SpecificEnergy[month]/(self.Deviations[month]*sqrt(2*pi)))*exp(-((TimeStep-self.NoonTime[month])**2)/(2*self.Deviations[month]**2))*erf(-self.Alpha[0]*(TimeStep-self.SunriseTime[month])/(self.Deviations[month]*sqrt(2)))*erf(-self.Alpha[1]*(TimeStep-self.SunsetTime[month])/(self.Deviations[month]*sqrt(2))) if logical_and(TimeStep>self.SunriseTime[month], TimeStep<self.SunsetTime[month]) else 0) for TimeStep in self.TimeStamp] for month in range(len(self.SpecificEnergy))]).transpose()

    def EstimateDCPowerOutput(self, DCPeakPower = 1):
        from numpy import append

        self.SunriseTime = self.Transform2DecimalHour(self.SunriseTime)
        self.SunsetTime = self.Transform2DecimalHour(self.SunsetTime)
        self.NoonTime = self.Transform2DecimalHour(self.NoonTime)
        self.PVEnvelopes = self.GaussPVEnvelopes()*DCPeakPower
        self.Power = array([])
        
        for month in range(12):
            if month == 0:
                for day in range(31):
                    self.Power = append(self.Power, self.PVEnvelopes[:,0])
                
            elif month == 1:
                for day in range(28):
                    self.Power = append(self.Power, self.PVEnvelopes[:,1])
                
            elif month == 2:
                for day in range(31):
                    self.Power = append(self.Power, self.PVEnvelopes[:,2])
                
            elif month == 3:
                for day in range(30):
                    self.Power = append(self.Power, self.PVEnvelopes[:,3])
                
            elif month == 4:
                for day in range(31):
                    self.Power = append(self.Power, self.PVEnvelopes[:,4])
                
            elif month == 5:
                for day in range(30):
                    self.Power = append(self.Power, self.PVEnvelopes[:,5])
                
            elif month == 6:
                for day in range(31):
                    self.Power = append(self.Power, self.PVEnvelopes[:,6])
                
            elif month == 7:
                for day in range(31):
                    self.Power = append(self.Power, self.PVEnvelopes[:,7])
                
            elif month == 8:
                for day in range(30):
                    self.Power = append(self.Power, self.PVEnvelopes[:,8])
                
            elif month == 9:
                for day in range(31):
                    self.Power = append(self.Power, self.PVEnvelopes[:,9])
                
            elif month == 10:
                for day in range(30):
                    self.Power = append(self.Power, self.PVEnvelopes[:,10])
                
            else:
                for day in range(31):
                    self.Power = append(self.Power, self.PVEnvelopes[:,11])
      

def BESS_perm_min(SoC, Capacity_BESS, SoCmax, P_BESS_max, dt = 0.25):
    from numpy import clip
    
    return clip(Capacity_BESS*(SoC - SoCmax)/dt, -P_BESS_max, P_BESS_max)
    
    
def BESS_perm_max(SoC, Capacity_BESS, SoCmin, P_BESS_max, dt = 0.25):
    from numpy import clip
    
    return clip(Capacity_BESS*(SoC - SoCmin)/dt, -P_BESS_max, P_BESS_max)

class BESS:
    def __init__(self, Capacity = 0, MaxPower = 0, Char_eff = 0.92, Dischar_eff = 0.92, Self_discharge = 0, MinSoC = 0.1, MaxSoC = 0.9, Power = [], Energy = [], SoC = []):
        self.Capacity = Capacity                    # In kWh
        self.MaxPower = MaxPower                    # In kW
        self.Char_eff = 1/Char_eff                  # In %
        self.Dischar_eff = 1/Dischar_eff            # In %
        self.Self_discharge = (1-Self_discharge)    # In %
        self.MinSoC = MinSoC                        # In %
        self.MaxSoC = MaxSoC                        # In %
        self.Power =  Power                         # In kW
        self.Energy = Energy                        # In kWh
        self.SoC = SoC                              # In %

class PeakShaving:
    def __init__(self, PV, BESS, Load, MaxPower = 0, StartTime = 0, EndTime = 0, TimeStep = 0.25):
        self.PV = PV
        self.BESS = BESS
        self.Load = Load
        self.Grid = []
        self.MaxPower = MaxPower
        self.StartTime = StartTime
        self.EndTime = EndTime
        self.TimeStep = TimeStep
    
    def Initialize(self, InitialSoC = 0, InitialPower = 0):
        
        self.BESS.Power = self.PV*0 
        self.BESS.Power[0] = InitialPower
        self.BESS.Energy =  self.PV*0     
        self.BESS.Energy[0] = InitialSoC*self.BESS.Capacity
        self.BESS.SoC =  self.PV*0 
        self.BESS.SoC[0] = InitialSoC
        self.Grid =  self.PV*0 
        self.Grid[0] = self.Load[0] - self.PV[0] - InitialPower 

        
    def Simulate(self):
        from math import floor    
        
        # Consider the power of the BESS inverter.
        
        for t_step in range(floor((self.EndTime-self.StartTime)/self.TimeStep)-1):
            step = t_step+1
            
            if self.BESS.SoC[step-1] >= self.BESS.MaxSoC:                  # Battery can only discharge

                if (self.Load[step]-self.PV[step]) <= self.MaxPower:        # No peakshaving needed
                    self.BESS.Power[step] = 0
                    self.BESS.Energy[step] = self.BESS.Energy[step-1]*self.BESS.Self_discharge
                    self.BESS.SoC[step] = self.BESS.SoC[step-1]
                    self.Grid[step] = self.Load[step] - self.PV[step] 

                
                else:                                                       # Peakshaving needed
                    if self.Load[step] - self.PV[step] - self.MaxPower  <= BESS_perm_max(self.BESS.SoC[step-1], Capacity_BESS = self.BESS.Capacity, SoCmin = self.BESS.MinSoC, P_BESS_max = self.BESS.MaxPower):    # Below the BESS max power
                        self.BESS.Power[step] = self.Load[step] - self.PV[step] - self.MaxPower
                        self.BESS.Energy[step] = self.BESS.Energy[step-1] - self.BESS.Power[step]*self.TimeStep*self.BESS.Dischar_eff
                        self.BESS.SoC[step] = self.BESS.Energy[step] /  self.BESS.Capacity
                        self.Grid[step] = self.MaxPower
                        
                    else:                                                       # Above the BESS max power
                        self.BESS.Power[step] = BESS_perm_max(self.BESS.SoC[step-1], Capacity_BESS = self.BESS.Capacity, SoCmin = self.BESS.MinSoC, P_BESS_max = self.BESS.MaxPower)
                        self.BESS.Energy[step] = self.BESS.Energy[step-1] - self.BESS.Power[step]*self.TimeStep*self.BESS.Dischar_eff
                        self.BESS.SoC[step] = self.BESS.Energy[step] /  self.BESS.Capacity
                        self.Grid[step] = -self.BESS.MaxPower - self.PV[step] + self.Load[step]              
            
            elif self.BESS.SoC[step-1] > self.BESS.MinSoC:                  # Battery can charge and discharge
                
                if (self.Load[step]-self.PV[step]) <= self.MaxPower:        # No peakshaving needed
                    
                    if self.Load[step] >= self.PV[step]:                    # PV below demanded load
                        self.BESS.Power[step] = 0
                        self.BESS.Energy[step] = self.BESS.Energy[step-1]*self.BESS.Self_discharge
                        self.BESS.SoC[step] = self.BESS.SoC[step-1]
                        self.Grid[step] = self.Load[step] - self.PV[step]                      

                    else:                                                   # Surplus of PV power
                        if self.PV[step] - self.Load[step] <= -BESS_perm_min(self.BESS.SoC[step-1], Capacity_BESS = self.BESS.Capacity, SoCmax = self.BESS.MaxSoC, P_BESS_max = self.BESS.MaxPower):    # Below the BESS max power
                            self.BESS.Power[step] = self.Load[step] - self.PV[step]
                            self.BESS.Energy[step] = self.BESS.Energy[step-1] - self.BESS.Power[step]*self.TimeStep*self.BESS.Char_eff
                            self.BESS.SoC[step] = self.BESS.Energy[step] /  self.BESS.Capacity
                            self.Grid[step] = 0

                            
                        else:                                                       # Above the BESS max power
                            self.BESS.Power[step] = BESS_perm_min(self.BESS.SoC[step-1], Capacity_BESS = self.BESS.Capacity, SoCmax = self.BESS.MaxSoC, P_BESS_max = self.BESS.MaxPower)
                            self.BESS.Energy[step] = self.BESS.Energy[step-1] - self.BESS.Power[step]*self.TimeStep*self.BESS.Char_eff
                            self.BESS.SoC[step] = self.BESS.Energy[step] /  self.BESS.Capacity
                            self.Grid[step] = self.BESS.MaxPower - self.PV[step] + self.Load[step]

                            
                else:                                                       # Peakshaving needed
                    if self.Load[step] - self.PV[step] - self.MaxPower <= BESS_perm_max(self.BESS.SoC[step-1], Capacity_BESS = self.BESS.Capacity, SoCmin = self.BESS.MinSoC, P_BESS_max = self.BESS.MaxPower):    # Below the BESS max power
                        self.BESS.Power[step] = self.Load[step] - self.PV[step] - self.MaxPower
                        self.BESS.Energy[step] = self.BESS.Energy[step-1] - self.BESS.Power[step]*self.TimeStep*self.BESS.Dischar_eff
                        self.BESS.SoC[step] = self.BESS.Energy[step] /  self.BESS.Capacity
                        self.Grid[step] = self.MaxPower

                        
                    else:                                                       # Above the BESS max power
                        self.BESS.Power[step] = BESS_perm_max(self.BESS.SoC[step-1], Capacity_BESS = self.BESS.Capacity, SoCmin = self.BESS.MinSoC, P_BESS_max = self.BESS.MaxPower)
                        self.BESS.Energy[step] = self.BESS.Energy[step-1] - self.BESS.Power[step]*self.TimeStep*self.BESS.Dischar_eff
                        self.BESS.SoC[step] = self.BESS.Energy[step] /  self.BESS.Capacity
                        self.Grid[step] = -self.BESS.MaxPower - self.PV[step] + self.Load[step]  

            
            else: # self.BESS.SoC[1,i] <= self.BESS.MinSoC:                 # Battery can only charge     
                
                if self.Load[step] >= self.PV[step]:                        # PV below demanded load
                    self.BESS.Power[step] = 0
                    self.BESS.Energy[step] = self.BESS.Energy[step-1]*self.BESS.Self_discharge
                    self.BESS.SoC[step] = self.BESS.SoC[step-1]
                    self.Grid[step] = self.Load[step] - self.PV[step] 
                     
                else:                                                       # Surplus of PV power
                    
                    if self.PV[step] - self.Load[step] <= -BESS_perm_min(self.BESS.SoC[step-1], Capacity_BESS = self.BESS.Capacity, SoCmax = self.BESS.MaxSoC, P_BESS_max = self.BESS.MaxPower):    # Below the BESS max power
                        self.BESS.Power[step] = self.Load[step] - self.PV[step]
                        self.BESS.Energy[step] = self.BESS.Energy[step-1] - self.BESS.Power[step]*self.TimeStep*self.BESS.Char_eff
                        self.BESS.SoC[step] = self.BESS.Energy[step] /  self.BESS.Capacity
                        self.Grid[step] = 0

                        
                        
                    else:                                                       # Above the BESS max power
                        self.BESS.Power[step] = BESS_perm_min(self.BESS.SoC[step-1], Capacity_BESS = self.BESS.Capacity, SoCmax = self.BESS.MaxSoC, P_BESS_max = self.BESS.MaxPower)
                        self.BESS.Energy[step] = self.BESS.Energy[step-1] - self.BESS.Power[step]*self.TimeStep*self.BESS.Char_eff
                        self.BESS.SoC[step] = self.BESS.Energy[step] /  self.BESS.Capacity
                        self.Grid[step] = self.BESS.MaxPower - self.PV[step] + self.Load[step]   

            
## As most of the data for the simulations is in .mat files, this function
## transforms a column of the .mat file into a list.

def mat2list(C):
   L1 = []
   for i in range(len(C)):
       j = C[i]
       L1.append(j[0])
   return L1

## As the equations use degrees instead of radians, these functins allow to 
## use degrees instead of radians.

def tand(x):
    return float(sp.tan(x * sp.pi / 180))

def sind(x):
    return float(sp.sin(x * sp.pi / 180))

def cosd(x):
    return float(sp.cos(x * sp.pi / 180))

def atand(x):
    return float((180/sp.pi)*sp.atan(x))

def asind(x):
    return float((180/sp.pi)*sp.asin(x))

def acosd(x):
    return float((180/sp.pi)*sp.acos(x))

## Julian date.

###############################################################################
############################   Classes definition  ############################

## This class storages the meteorological data, stored in a .mat file,
## containing atributes for DHI, DNI, GHI, ambient temperature and wind speed.
        
class met():
    
    def __init__(self):
        self.DHI = []
        self.DNI = []
        self.GHI = []
        self.T_amb = []
        self.u = []
        
    def load_data(self, doc=''):
        mat=scipy.io.loadmat(doc)
        self.DHI = array([[element for element in upperElement] for upperElement in mat['DHI']])
        self.DNI = array([[element for element in upperElement] for upperElement in mat['DNI']])
        self.GHI = array([[element for element in upperElement] for upperElement in mat['GHI']])
        self.T_amb = array([[element for element in upperElement] for upperElement in mat['T_amb']])
        self.u = array([[element for element in upperElement] for upperElement in mat['u']])

## This class storages the load data.

class load():
 
    def __init__(self):
        self.L = []
    
    def load_data(self, doc=''):
        mat=scipy.io.loadmat(doc)
        self.L = array([[element for element in upperElement] for upperElement in mat['L_CR_1h_5kWh']])

class Victor_data():
    def __init__(self):
        self.JDD = []
        self.P_pv = []

    def load_data(self, doc=''):
        import scipy.io
        
        mat=scipy.io.loadmat(doc)
        self.JDD = array([[element for element in upperElement] for upperElement in mat['JDD']])
        self.P_pv = array([[element for element in upperElement] for upperElement in mat['P_pv']])


###############################################################################
###### INFORMATION ABOUT THE MODULE INCLINATION AND AZIMUTH ######
## This values normally rely on previously identified optimum tilt and azimth
## obtained thorugh a annual energy production analyss.  

class PV_system():
    
    def __init__(self):
        
        # Module orientation parameters
        self.Opt_theta = 0
        self.Opt_Am = 0
        self.alpha = 0
        
        # Module datasheet parameters
        self.eff_datasheet = 0 #given by the datasheet
        self.eff = 0
        self.C_eff = 0 # thermal coeficient 
        self.thao_alpha = 0
        self.Tm_NOCT = 0 #Temperature of the PV module at NOCT, data sheet  
        self.Ta_NOCT = 0
        self.Voc_STC = 0
        self.Isc_STC = 0
        self.Vmpp_STC = 0
        self.Impp_STC = 0
        self.Pmmp_STC = 0
        self.FF_STC = 0
        self.Area = 0
        
        self.Gm_STC = 1000
        self.Gm_NOCT = 800
        self.nn = 1.5
        self.k_b = 1.36e-23
        self.e = 1.602e-19
        
        self.T_STC = 25
        self.T_STC_K = 273.15+self.T_STC
        
        # Inverter/charger datasheet parameters
        self.eff_mppt = 0
        self.eff_inverter = 0
        
        # Battery datasheet parameters
        self.V_nom = 0
        self.n_serie = 0
        self.n_parallel = 0
        self.eff_battery = 0 # battery efficiency, considered constant, CAN BE IMPROVED
        self.SoC_lowerlimit = 0 # lower limit to define the range of operation of the battery
        self.SoC_upperlimit = 0 # upper limit to define the range of operation of the battery   
        self.eff_CC = 0       
        self.V_battery = 0 #total pack voltage
        self.cell_cap = 0
        self.Cap_battery = 0 # battery nominal capacity in Ah is the capacity of the pouch cell selected from A123
        self.E_batt_cap = 0 # nominal battery energy 
        self.lowerlimit_batt = 0 # lower allowed battery enegy
        self.upperlimit_batt = 0 # maximum allowed battery enegy
        self.usable = 0  
        
    def load_module_orientation(self, Opt_theta=10.5, Opt_Am=200, alpha=0.2):  # The default parameters are based on Costa Rican conditions.
        self.Opt_theta = Opt_theta
        self.Opt_Am = Opt_Am
        self.alpha = alpha
        
        
    def load_module_datasheet(self, eff_datasheet, C_eff, thao_alpha, Tm_NOCT, Ta_NOCT, Voc_STC, Isc_STC, Vmpp_STC, Impp_STC, Area):
        
        self.eff_datasheet = eff_datasheet
        self.eff = self.eff_datasheet
        self.C_eff = C_eff
        self.thao_alpha = thao_alpha
        self.Tm_NOCT = Tm_NOCT   
        self.Ta_NOCT = Ta_NOCT
        self.Voc_STC = Voc_STC
        self.Isc_STC = Isc_STC
        self.Vmpp_STC = Vmpp_STC
        self.Impp_STC = Impp_STC
        self.Area = Area
             
        self.Pmmp_STC = self.Impp_STC*self.Vmpp_STC
        self.FF_STC = self.Pmmp_STC/(self.Voc_STC*self.Isc_STC)
    
    def load_inverter_datasheet(self, eff_mppt, eff_inverter):
        self.eff_mppt = eff_mppt # efficiency of the mppt process
        self.eff_inverter = eff_inverter # efficiency of the inverter, it can be claculated dynamically 
        
    def load_battery_datasheet(self, V_nom, n_serie, cell_cap, n_parallel, eff_battery, SoC_lowerlimit, SoC_upperlimit, eff_CC=0):
        self.V_nom = V_nom
        self.n_serie = n_serie
        self.n_parallel = n_parallel
        self.cell_cap = cell_cap
        self.eff_battery = eff_battery # battery efficiency, considered constant, CAN BE IMPROVED
        self.SoC_lowerlimit = SoC_lowerlimit # lower limit to define the range of operation of the battery
        self.SoC_upperlimit = SoC_upperlimit # upper limit to define the range of operation of the battery
        
        self.eff_CC = eff_CC
        
        self.V_battery = self.V_nom*self.n_serie #total pack voltage
        self.Cap_battery = self.cell_cap*self.n_parallel # battery nominal capacity in Ah is the capacity of the pouch cell selected from A123
        self.E_batt_cap = self.Cap_battery*self.V_battery # nominal battery energy 
        self.lowerlimit_batt = self.SoC_lowerlimit*self.E_batt_cap # lower allowed battery enegy
        self.upperlimit_batt = self.SoC_upperlimit*self.E_batt_cap # maximum allowed battery enegy
        self.usable = (self.SoC_upperlimit-self.SoC_lowerlimit)*self.E_batt_cap
        
###############################################################################
###### MODEL ######
        
class NREL_PV_model():
    
    def __init__(self):
        
        # Geographical data
        self.L = 0
        self.l = 0
        self.N_m = 0
        self.Time_Zone = ''
        
        # System variables
        
        self.Tm_DB_it = 0 
        self.Isc_Gm = 0
        self.Voc_Gm = 0
        self.Pmmp_Gm = 0
        self.eff_Gm = 0
        self.register_eff = 0
        
        # Time variables
        self.s = 0
        self.e = 0
        self.dt = 0    
        self.ht1 = 0
        self.ht2 = 0
        
        # Simulation start point
        self.year_s = 0
        self.month_s = 0
        self.day_s = 0
        self.hour_s = 0
        self.minute_s = 0
        self.second_s = 0
        
        # Simulation end point
        self.year_e = 0
        self.month_e = 0
        self.day_e = 0
        self.hour_e = 0
        self.minute_e = 0
        self.second_e = 0
  
        
        # Inicialization of matrixes
        self.P_pv = 0
        self.P_pv_syst = array([]) #array([0]*365*24*4)
        self.cos_AOI = 0
        self.Gdirect = 0
        self.Gdiffuse = 0
        self.Galbedo = 0
        self.Gm = 0



    def sim_init(self, met_data, PV_system, step, year_s, month_s, day_s, hour_s, minute_s, second_s, year_e, month_e, day_e, hour_e, minute_e, second_e, L, l, N_m, Time_zone):        
        from pandas import Timestamp
        from math import floor, log
        from datetime import datetime, timedelta
        from numpy import append
        
        
        # Geographical data
        self.L = L
        self.l = l
        self.N_m = N_m
        self.Time_Zone = Time_zone
        
        # Time initialiation
        
        self.step = step
      
               
        # Simulation range
        # normally a complete year is used (365 days) to include the effect of the seasons in the calculations
    
        #Se define el formato de la fecha 
        #Es necesario especificar la zona horaria correspondiente, por ejemplo: para los datos solares de Delft, se elige
        #'Europe/Amsterdam'. 
        #Para conocer las zonas horarias existentes se utiliza el comando T = timezones('X'), donde la letra X 
        #se puede reemplazar por alguna de las siguientes localidades: Africa, America, Antarctica, Arctic, Asia, 
        #Atlantic, Australia, Etc, Europe, Indian, Pacific.
        
        self.s = datetime(year_s, month_s, day_s, hour_s, minute_s, second_s) #formating the staing date
        self.e = datetime(year_e, month_e, day_e, hour_e, minute_e, second_e) #formating the ending date        
        self.dt = self.e-self.s
        self.time_delta = timedelta(hours=1)
        
        #Tiempo de inicio de la simulacion expresado en horas
        self.ht1 = hour_s+minute_s/60+second_s/3600
        self.ht2 = self.ht1+self.dt.days*24
    
    
        #for i in range(int(self.ht2)):
        #self.t_per_step = [0]*365*24
        for i in range(365*24):
            ts = Timestamp(year = self.s.year, month = self.s.month, day = self.s.day, hour = self.s.hour, minute = self.s.minute, second = self.s.second, tz = Time_zone)
            D=ts.to_julian_date()-2451545+0.25
            
            
            # 2. Obtaining the mean longitude q and the mean anomaly g
            # note: q and g must be normalised
            q_i=(280.459+0.98564736*D)/360 #to know how many complete turns 
            q_n=(q_i-floor(q_i))*360 #gives q_normilised
            

            g_i=(357.529+0.98560028*D)/360; #to know how many complete turns 
            g_n=(g_i-floor(q_i))*360 #gives q_normilised
            
            #3. Calculating the eliptic longitude lamda_s
            lamda_s=q_n+1.915*sind(g_n)+0.02*sind(2*g_n)
            
            #4. Calculating axial tilt (empsilon) 
            empsilon=23.429-0.00000036*D
            
            #5. Getting Greenwich mean sidereal time (GMST)
            T=D/36525 # normally not taking into account
            #GMST=18.697374558+24.06570982441908*D+0.000026*T**2 #here D is given in days
            #Note: here, the term 0.000026*T^2 is insignificant 
            #GMST_i=GMST/24 #to normalise by 24h  
            GMST_i=(18.697374558+24.06570982441908*D+0.000026*T**2)/24 #to normalise by 24h 
            GMST_n=(GMST_i-floor(GMST_i))*24 #gives GMST_normalised
            
            #6. Estimating local mean sidearal time (LMST or theta_L)
            LMST=GMST_n*15+self.l
            
            #7. Obtaining azimuth (A_s) and altitud (a_s) of the sun 
            
            N_tan_A_s=(-sind(LMST)*cosd(lamda_s)+cosd(LMST)*cosd(empsilon)*sind(lamda_s)) # numerator for calculating tangent of the azimuth
            D_tan_A_s=(-sind(self.L)*cosd(LMST)*cosd(lamda_s)-(sind(self.L)*sind(LMST)*cosd(empsilon)-cosd(self.L)*sind(empsilon))*sind(lamda_s)) # denominator for calculating tangent of the azimuth
            tan_A_s=N_tan_A_s/D_tan_A_s     # here, the tangent of azimuth is calculated
            
            sin_a_s=cosd(self.L)*cosd(LMST)*cosd(lamda_s)+(cosd(self.L)*sind(LMST)*cosd(empsilon)+sind(self.L)*sind(empsilon))*sind(lamda_s) # here, the sin of altitude is calculated
            
            a_s=asind(sin_a_s) # Finally, this is the altitude of the sun
             
            #7.1 Getting the azimuth of the sun is not and straight foward task. To know its value, a couple of onditions must be evaluated first  
            if D_tan_A_s > 0 and N_tan_A_s > 0:
                A_s=atand(tan_A_s) # Finally, this is the azimuth of the sun if the numerator is possitive as well as the denominator
            elif D_tan_A_s < 0: 
                A_s=atand(tan_A_s) + 180 # Finally, this is the azimuth of the sun if the denominator is negative
            elif D_tan_A_s > 0 and N_tan_A_s < 0:
                A_s=atand(tan_A_s) + 360 # Finally, this is the azimuth of the sun if the denominator is positive but the numerator negative
            
            
            if a_s < 0:
                a_s = 0
                       
            
            #Se actualiza la fecha para la siguiente iteracion
            self.s += self.time_delta
            
                        
            # Direct Irradiance
            
            self.cos_AOI = cosd(90-PV_system.Opt_theta)*cosd(a_s)*cosd(PV_system.Opt_Am-A_s)+sind(90-PV_system.Opt_theta)*sind(a_s) # Calculate the cosine of the angle of incidence for every hour of the year for a specific module tilt and orientation            
            self.Gdirect = met_data.DNI[i]*self.cos_AOI #To calculate the direct irradiance 
            if self.Gdirect < 0:   # This makes al the negative Direct Irradiance on the PV module as zero.
                self.Gdirect = 0


            
            # Isotropic Diffuse Irradiance
            SVF = (1+cosd(PV_system.Opt_theta))/2 #to calculate the sky view factor
            self.Gdiffuse = met_data.DHI[i]*SVF
            
            # Irradiance due to effect of albedo
            self.Galbedo = met_data.GHI[i]*PV_system.alpha*(1-SVF)            
                
            # Total Irradinace
            self.Gm = self.Gdirect+self.Gdiffuse+self.Galbedo
            
 
            if self.Gm>1: # considering the effect of irradiance               
                Isc_Gm=PV_system.Isc_STC*self.Gm/PV_system.Gm_STC # just considering Gm for getting short circuit current
                Voc_Gm=PV_system.Voc_STC+PV_system.nn*PV_system.k_b*PV_system.T_STC_K/PV_system.e*log(self.Gm/PV_system.Gm_STC) # just considering Gm for getting open circuit voltage           
                Pmmp_Gm=PV_system.FF_STC*Isc_Gm*Voc_Gm           
                eff_Gm=Pmmp_Gm/(self.Gm*PV_system.Area) # efficiency just including effect of irradiance 

    
                for j in range (5):
                    Tm_DB_it=met_data.T_amb[i]+(PV_system.Tm_NOCT-PV_system.Ta_NOCT)*(self.Gm/PV_system.Gm_NOCT)*(9.5/(5.7+3.8*met_data.u[i]))*(1-(PV_system.eff/(PV_system.thao_alpha))) # using Duffie - Beckmann model
                    eff=eff_Gm*(1+PV_system.C_eff*(Tm_DB_it-(PV_system.T_STC))) # this finally take into account the thermal and irradiation effect over the efficiency of the pv module
    
            elif self.Gm<=1:  # 1 W/m^2                 
                Tm_DB_it=met_data.T_amb[i] # assign ambient temperature when there is no irradiance or it is minimal
                eff=0;    

            self.P_pv = self.Gm*PV_system.Area*eff # power produced per pv module            
            for j in range(4):        
                self.P_pv_syst = append(self.P_pv_syst, self.P_pv*self.N_m*PV_system.eff_mppt/1000)




###############################################################################
###### Costs ######

def year2months(yearly_data, dt, leapyear = False):
    dt = int(24/dt)
    
    if leapyear == False:
        months = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]                   
    else:
        months = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
    
    for i in range(12):
        monthly_data = [yearly_data[int(months[i]*dt):int(months[i+1]*dt)] for i in range(12)]
        
    return monthly_data 
    

def yearlycosts(yearly_data, timeframes, hour_costs, dt = 0.25):
    costs_list = []
    hour = 0
    timeframe_index = 0
    
    for time_step in yearly_data:
        
        if hour < timeframes[timeframe_index]:
            costs_list.append(hour_costs[timeframe_index])
            
        else:
            if hour == 24:
                hour = 0
                timeframe_index = 0
            else:
                timeframe_index +=1
            costs_list.append(hour_costs[timeframe_index])
                                    
        hour += dt
        
    return costs_list
         
###############################################################################