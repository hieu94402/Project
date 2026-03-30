import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
start = time.time()
# define constant
v_f = 1/60 # feeding speed, min/s
v_d = 100/60 # drawing speed, min/s
# L # furnace length, m, defined later
L_T_max = 0.1 # position at which temp peaks
# T_max # max temperature, degree
lambda_0 = 200 # initial periodicity, μm
alpha = 14300 # temp profile parameter formula 0

# case0:
L = 0.2
T_max = np.linspace(150, 200, 30) # define max temperature range
def _T(z, T_max): # temperature profile
    return T_max - alpha * (z-L_T_max) * (z-L_T_max)
def _eta(z, T_max): # viscosity profile
    return np.exp(22493 / (_T(z, T_max) + 273.15) - 35.287)
def _v(z, T_max): # velocity profile
    _v_denom, _ = quad(lambda g: 1.0 / _eta(g, T_max), 0, L)
    _v_num, _ = quad(lambda g: 1.0 / _eta(g, T_max), 0, z)
    return np.exp(np.log(v_f) + (_v_num/_v_denom) * np.log(v_d/v_f))
def _gamma(z, T_max): # surface tension profile
    return 49.2 - 0.06 * (_T(z, T_max) - 20)
def _lambda(z, T_max): # periodicity profile
    return lambda_0 * np.sqrt(v_f / _v(z,T_max))
def _tau(z, T_max): # characteristic time
    return _eta(z, T_max) * 1e-6 * _lambda(z, T_max) / (3.14 * _gamma(z, T_max))
def _fsh(T_max): #shape factor
    _fsh_func, _ = quad(lambda g: -1.0 / (_tau(g, T_max) * _v(g, T_max)), 0, L)
    return np.exp(_fsh_func)
# modify here
y0 = np.vectorize( _fsh )
# _proc_values = _proc_vec( T_max )
# plt.plot(T_max, _proc_values, label="case0: parabola-like decay", color='C0')
cols = []
cols.append(T_max)
cols.append(y0(T_max))

# case1:
L = 0.3 # re-define furnace length
exp_diff = 46
def _T(z, T_max): 
    y_piecewise = np.where(
        z <= L_T_max,
        T_max - alpha * (z-L_T_max) * (z-L_T_max),
        T_max * np.exp(- exp_diff * (z-L_T_max)**2))
    return y_piecewise
y1 = np.vectorize( _fsh )
cols.append(y1(T_max))

# case2:
exp_diff = 35
y2 = np.vectorize( _fsh )
cols.append(y2(T_max))

# case3:
exp_diff = 100
y3 =np.vectorize( _fsh )
cols.append(y3(T_max))

# plt.legend()
# plt.show() # show figure

data = np.column_stack(cols)
np.savetxt('C:\\Users\\hieu9\\OneDrive\\Máy tính\\One\\[Project] Shape preservation and stress relaxation\\python calculation\\data\\results-fsh.csv', data, delimiter=',',)

end = time.time()
runtime = end - start
# print(f"Runtime: {runtime:.3f} s")