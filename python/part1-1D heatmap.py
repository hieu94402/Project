import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime
from scipy.integrate import simpson

# define constant
v_f = 1/6000 # feeding speed, mm/s
v_d = 1000/6 # drawing speed, mm/s
# L1 # begin coordinate, m
# L2 # end coordinate, m
# position at which temp peaks is z = 0
# T_max # max temperature, degree
lambda_0 = 200 # initial periodicity, μm
alpha = 14300 # temp profile parameter formula 0

# def _T(z): # temperature profile
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

T_max = 170
L1 = -0.1
L2 = 0.5
z = np.linspace(L1, L2, 1000) # define horizontal coor.
chi = 7
def _T(z): # temperature profile
    y_piecewise = np.where(
        z < 0, 
        T_max - alpha * z **2,
        T_max * np.exp(- chi * z ** 2))
    return y_piecewise
y = np.vectorize( _T )
T_values = _T(z)
# reshape to (1, N) for imshow — a single-row 2D array
T_strip = T_values.reshape(1, -1)

fig, axes = plt.subplots(
    nrows=2, ncols=1,
    figsize=(10, 3),
    gridspec_kw={'height_ratios': [1, 0.15]},  # main heatmap + thin colorbar row
    constrained_layout=True
)
# --- main heatmap strip ---
im = axes[0].imshow(
    T_strip,
    aspect='auto',
    cmap='RdYlBu_r',          # change to 'plasma', 'hot', 'RdYlBu_r', etc.
    extent=[L1, L2, 0, 1],  # [x_min, x_max, y_min, y_max]
    origin='lower'
)

# overlay a line showing exact T profile for reference (optional)
ax2 = axes[0].twinx()
ax2.plot(z, T_values, color='white', linewidth=1, linestyle='--', alpha=0.6, label='T(z)')
ax2.set_ylabel("Temperature (°C)", color='white')
ax2.tick_params(axis='y', colors='white')
ax2.yaxis.label.set_color('white')

axes[0].set_xlabel("z coordinate (m)")
axes[0].set_yticks([])          # hide y-axis ticks on heatmap
axes[0].set_title("Temperature Profile — Heatmap")

# --- colorbar below ---
cbar = fig.colorbar(im, cax=axes[1], orientation='horizontal')
cbar.set_label("Temperature (°C)")

plt.show()
