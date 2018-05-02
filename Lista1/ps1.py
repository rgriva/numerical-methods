# Metodos Numericos - EPGE/FGV 2018
# Instructor: Cezar Santos
# Problem Set 1 - Raul Guarini Riva

import numpy as np      # Linear Algebra Library
import time         # Just to compute execution time
from scipy.stats import norm
import gc

start_time = time.time()
print('Starting Problem Set 1...')
# Defining the Tauchen Discretization function

# Something VERY IMPORTANT: Python indexing starts from zero and ends on N-1 for a size-N array


def tauchen_ar1(mu, rho, sigma2, N, m):
    sigma = np.sqrt(sigma2)
    theta_N = m * sigma / np.sqrt((1 - rho**2))
    theta_1 = - m * sigma / ((1 - rho**2)**0.5)
    delta = (theta_N - theta_1) / (N - 1)

    grid = np.linspace(start=theta_1, stop=theta_N, num=N) + mu
    J, I = np.meshgrid(grid, grid)
    u_bound = (1 / sigma) * (J + (delta / 2) - (1 - rho) * mu - rho * I)
    l_bound = (1 / sigma) * (J - (delta / 2) - (1 - rho) * mu - rho * I)
    P = norm.cdf(u_bound) - norm.cdf(l_bound)

    # Fixing terminal points
    P[:, 0] = norm.cdf((1 / sigma) * (grid[0] - (1 - rho) * mu - rho * I[:, 0] + delta / 2))
    P[:, N - 1] = 1 - norm.cdf((1 / sigma) * (grid[N - 1] - (1 - rho)
                                              * mu - rho * I[:, 0] - delta / 2))
    return grid, P

# Defining the Rouwenhorst Discretization function


def rou(mu, rho, sigma2, N):
    theta_N = np.sqrt(sigma2 / (1 - rho**2)) * np.sqrt(N - 1)
    theta_1 = - theta_N
    grid = np.linspace(start=theta_1, stop=theta_N, num=N) + mu
    p = (1 + rho) / 2
    P_init = np.array([[p, 1 - p], [1 - p, 1]])
    for n in range(3, N + 1):
        stack1 = np.vstack((np.hstack((P_init, np.zeros((n - 1, 1)))), np.zeros((1, n))))
        stack2 = np.vstack((np.hstack((np.zeros((n - 1, 1)), P_init)), np.zeros((1, n))))
        stack3 = np.vstack((np.zeros((1, n)), np.hstack((P_init, np.zeros((n - 1, 1))))))
        stack4 = np.vstack((np.zeros((1, n)), np.hstack((np.zeros((n - 1, 1)), P_init))))
        P_new = p * stack1 + (1 - p) * stack2 + (1 - p) * stack3 + p * stack4
        # Updating the loop
        P_init = P_new

    # Normalizing conditional measures
    P = P_new / P_new.sum(axis=1, keepdims=True)
    return grid, P


# Item 1
# Calibration
mu = 0
rho = 0.95
sigma = 0.007
sigma2 = sigma**2
N = 9
m = 3      # Scaling parameter

grid_tauchen, P_tauchen = tauchen_ar1(mu, rho, sigma2, N, m)
print('Tauchen discretization done.')

# Item 2
grid_rou, P_rou = rou(mu, rho, sigma2, N)
print('Rouwenhorst discretization done.')

# Item 3
burn = 100
periods = 10000

# Drawing the shock sequence
shocks = norm.rvs(loc=0, scale=sigma, size=(periods + burn, 1))

# Simulating the continuous path
cont_path = np.zeros((periods + burn, 1))
cont_path[0] = mu

for t in range(1, periods + burn):
    cont_path[t] = (1 - rho) * mu + rho * cont_path[t - 1] + shocks[t]
# Deleting the burn periods, so initial condition is stochastic
cont_path = cont_path[burn - 1:-1]
print('Continuous Path simulated.')

# Mapping the normal shocks on uniform ones
uni_shocks = norm.cdf(shocks, loc=0, scale=sigma)

# Simulating Tauchen path
tauchen_path = np.zeros((periods + burn, 1))
state = np.floor(N / 2).astype('uint8')
tauchen_path[0] = grid_tauchen[state]
for t in range(1, periods + burn):
    cum_sum = np.cumsum(P_tauchen[state, :])
    ind = sum(uni_shocks[t] >= cum_sum)
    state = ind
    tauchen_path[t] = grid_tauchen[state]

tauchen_path = tauchen_path[burn - 1:-1]
print('Tauchen Path simulated.')

# Simulating Rouwenhorst path
rou_path = np.zeros((periods + burn, 1))
state = np.floor(N / 2).astype('uint8')
rou_path[0] = grid_tauchen[state]
for t in range(1, periods + burn):
    cum_sum = np.cumsum(P_rou[state, :])
    ind = sum(uni_shocks[t] >= cum_sum)
    state = ind
    rou_path[t] = grid_rou[state]

rou_path = rou_path[burn - 1:-1]
print('Rouwenhorst Path simulated.')
gc.collect()

# Computing MSE
MSE_tauchen = sum((tauchen_path - cont_path)**2) / periods
MSE_rou = sum((rou_path - cont_path)**2) / periods
print('Mean Square Error of Tauchen Discretization: {:f}'.format(float(MSE_tauchen)))
print('Mean Square Error of Rouwenhorst Discretization: {:f}'.format(float(MSE_rou)))

# Item 4 - Estimation
from sklearn.linear_model import LinearRegression   # This is a library to estimating linear models

fit_tauchen = LinearRegression().fit(X=tauchen_path[0:-2], y=tauchen_path[1:-1])
fit_rou = LinearRegression().fit(X=rou_path[0:-2], y=rou_path[1:-1])
print('True rho: {:0.3f}'.format(rho))
print('Estimated rho (Tauchen): {:f}'.format(float(fit_tauchen.coef_)))
print('Estimated rho (Rouwenhorst): {:f}'.format(float(fit_rou.coef_)))

print('Total execution time: {:.3f} seconds'.format(time.time() - start_time))

# Plotting stuff
import matplotlib.pyplot as plt
figure, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 6))
ax1.plot(range(1, periods + 1), cont_path, 'b', label='Continuous')
ax1.plot(range(1, periods + 1), tauchen_path, 'r', label='Tauchen')
ax1.legend()
ax1.grid(True)
ax1.set_title('Tauchen vs Continuous - N = {}, rho = {}'.format(N, rho))
ax1.set_xlabel('Time')

ax2.plot(range(1, periods + 1), cont_path, 'b', label='Continuous')
ax2.plot(range(1, periods + 1), tauchen_path, 'g-', label='Rouwenhorst')
ax2.grid(True)
ax2.set_title('Rouwenhorst vs Continuous - N = {}, rho = {}'.format(N, rho))
ax2.legend()
ax2.set_xlabel('Time')
plt.tight_layout()
plt.show()
