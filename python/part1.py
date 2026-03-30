import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
import matplotlib as mpl

# define constant
v_f = 1/60 # feeding speed, mm/s
v_d = 100/60 # drawing speed, mm/s
# L # furnace length, m, defined later
L_T_max = 0.1 # position at which temp peaks
T_max = 170 # max temperature, degree
lambda_0 = 200 # initial periodicity, μm
alpha = 14300 # temp profile parameter formula 0

# # preset plot
# plt.figure(figsize=(8,5))
# plt.grid(True)
# plt.xlabel("z coordinate (m)")
# mpl.rcParams['lines.marker'] = 'o'
# mpl.rcParams['lines.markersize'] = 1
# mpl.rcParams['lines.linestyle'] = '-'

# case0:
L = 0.2
z = np.linspace(0, L, 1000) # define horizontal coor. (furnace length, 0 eq to the top)
def _T(z): # temperature profile
    return T_max - alpha * (z-L_T_max) * (z-L_T_max)
def _eta(z): # viscosity profile
    return np.exp(22493 / (_T(z) + 273.15) - 35.287)
def _log_eta(z):
    return np.log(_eta(z))
def _v(z): # velocity profile
    _v_denom, _ = quad(lambda g: 1.0 / _eta(g), 0, L)
    _v_num, _ = quad(lambda g: 1.0 / _eta(g), 0, z)
    return np.exp(np.log(v_f) + (_v_num/_v_denom) * np.log(v_d/v_f))
def _gamma(z): # surface tension profile
    return 49.2 - 0.06 * (_T(z) - 20)
def _lambda(z): # periodicity profile (in micro-m)
    return lambda_0 * np.sqrt(v_f / _v(z))
def _tau(z): # characteristic time (log)
    return _eta(z) * 1e-6 * _lambda(z) / (3.14 * _gamma(z))
def _log_tau(z):
    return np.log(_tau(z))
# # modify here
y0 = np.vectorize( _v )
# y1_values = y1(z)
# plt.plot(z, y1_values, label="case0: parabola decay", color='C0')
cols = []
cols.append(z)
cols.append(y0(z))

# case1:
L = 0.3
z = np.linspace(0, L, 1000) # define horizontal coor. (furnace length, 0 eq to the top)
exp_diff = 46 # temp profile parameter formula 1
# the exponential parameter is set so that temperature at the end of the furnace approaches room-temp threshold
def _T(z): # temperature profile
    y_piecewise = np.where(
        z <= L_T_max,
        T_max - alpha * (z-L_T_max) * (z-L_T_max),
        T_max * np.exp(-exp_diff* (z-L_T_max)**2))
    return y_piecewise
# modify here
y1 = np.vectorize( _v )
cols.append(z)
cols.append(y1(z))

# case2:
exp_diff = 35
y2 = np.vectorize( _v )
cols.append(y2(z))

# case3:
exp_diff = 100
y3 = np.vectorize( _v )
cols.append(y3(z))

# plt.title("tau in log scale")
# plt.legend()
# plt.show() # show figure

# save to .csv file
data = np.column_stack(cols)
np.savetxt('C:\\Users\\hieu9\\OneDrive\\Máy tính\\One\\[Project] Shape preservation and stress relaxation\\python calculation\\data\\results-v.csv', data, delimiter=',')
