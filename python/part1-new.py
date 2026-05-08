import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime
from scipy.integrate import simpson

# define constant
v_f = 1/60 # feeding speed, mm/s
v_d = 100/60 # drawing speed, mm/s
# L1 # begin coordinate, m
# L2 # end coordinate, m
# position at which temp peaks is z = 0
# T_max # max temperature, degree
lambda_0 = 200 # initial periodicity, μm
alpha = 14300 # temp profile parameter formula 0

# preset plot
plt.figure(figsize=(8,5))
plt.grid(True)
plt.xlabel("z coordinate (m)")
mpl.rcParams['lines.marker'] = 'o'
mpl.rcParams['lines.markersize'] = 1
mpl.rcParams['lines.linestyle'] = '-'

# initiate .csv struc.
cols = []

# case0:
# def _T(z): # temperature profile
#     y_piecewise = np.where(
#         z > 0,
#         T_max - alpha * z **2,
#         T_max * np.exp(- chi * z ** 2))
#     return y_piecewise
def _eta(z): # viscosity profile
    return np.exp(22493 / (_T(z) + 273.15) - 35.287)
def _log_eta(z):
    return np.log(_eta(z))
def _v(z): # velocity profile
# using simpson method
    g_full = np.linspace(L1,L2,1000)
    g_part = np.linspace(L1,z,1000)
    denom = simpson(1.0 / _eta(g_full), x=g_full)
    numer = simpson(1.0 / _eta(g_part), x=g_part)
    return np.exp(np.log(v_f) + (numer/denom) * np.log(v_d/v_f))
# using quad integrator
    # _v_denom, _ = quad(lambda g: 1.0 / _eta(g), L1, L2)
    # _v_num, _ = quad(lambda g: 1.0 / _eta(g), L1, z)
    # return np.exp(np.log(v_f) + (_v_num/_v_denom) * np.log(v_d/v_f))
def _gamma(z): # surface tension profile
    return 49.2 - 0.06 * (_T(z) - 20)
def _lambda(z): # periodicity profile (in micro-m)
    return lambda_0 * np.sqrt(v_f / _v(z))
def _tau(z): # characteristic time (log)
    return _eta(z) * 1e-6 * _lambda(z) / (3.14 * _gamma(z))
def _log_tau(z):
    return np.log(_tau(z))
def shape_factor():
    integrand_vals = np.array([-1.0 / (_tau(g) * _v(g)) for g in z])
    fs = simpson(integrand_vals, x=z)
    return np.exp(fs)
# select ouput object
# case0:
T_max = 170
L1 = -0.5
L2 = 0.1
z = np.linspace(L1, L2, 1000) # define horizontal coor.
chi = 7
def _T(z): # temperature profile
    y_piecewise = np.where(
        z > 0, 
        T_max - alpha * z **2,
        T_max * np.exp(- chi * z ** 2))
    return y_piecewise
y0 = np.vectorize( _v )
y0_values = y0(z)
cols.append(z)
cols.append(y0(z))
# get shape factor value
fsh = shape_factor()
print(f"Shape factor fsh = {fsh:.6e}\n")
# set plot-figure
plt.plot(z, y0_values, color='C0', label='slow heating, fast quenching')

# case1:
L1 = -0.1
L2 = 0.5
z = np.linspace(L1, L2, 1000) # define horizontal coor.
def _T(z): # temperature profile
    y_piecewise = np.where(
        z < 0, 
        T_max - alpha * z **2,
        T_max * np.exp(- chi * z ** 2))
    return y_piecewise
y0 = np.vectorize( _v )
y0_values = y0(z)
cols.append(z)
cols.append(y0(z))
# get shape factor value
fsh = shape_factor()
print(f"Shape factor fsh = {fsh:.6e}\n")
# set plot-figure
plt.plot(z, y0_values, color='C1', label='fast heating, slow quenching')


# plt.title("unknown title")
# plt.legend()
plt.xlabel("z coordinate (m)")
plt.ylabel("unknown")
plt.legend()
plt.show() 
# # save to .csv file
# data = np.column_stack(cols)
# timestamp = datetime.now().strftime("%y%m%d-%H%M%S")
# filename = f'C:\\Users\\hieu9\\OneDrive\\Máy tính\\One\\[Project] Preservation & residual stress\\python calculation\\data\\excel\\~file~name'
# np.savetxt(filename, data, delimiter=',')
# print(f"Saved to: {filename}\n")