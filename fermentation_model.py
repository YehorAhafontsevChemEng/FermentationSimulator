import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

mu_max = float(input("Enter the maximum specific growth rate (mu_max): "))
Ks = float(input("Enter the half-saturation constant (Ks): "))
Yxs = float(input("Enter the yield coefficient (Yxs): "))
alpha = float(input("Enter the growth-associated product yield (alpha): "))
beta = float(input("Enter beta (g product per g cells per hour): "))
tau = float(input("Enter the lag time constant tau in hours (0 = no lag): "))

def mu(S):
    return mu_max * S / (Ks + S)

def lag_factor(t):
    if tau == 0:
        return 1.0
    return 1 - np.exp(-t / tau)

Yps = 0.5   # g product per g substrate used to make it

def model(t, y):
    X, S, P = y
    S = max(S, 0)
    sat = S / (Ks + S)                  # goes to 0 when substrate runs out
    growth = lag_factor(t) * mu(S)
    dXdt = growth * X
    dPdt = (alpha * growth + beta * sat) * X
    dSdt = -(1 / Yxs) * growth * X - (1 / Yps) * dPdt
    return [dXdt, dSdt, dPdt]

sol = solve_ivp(model, [0, 24], [0.1, 20, 0],
                t_eval=np.linspace(0, 24, 300),
                rtol=1e-8, atol=1e-10)

plt.plot(sol.t, sol.y[0], label="Biomass X")
plt.plot(sol.t, sol.y[1], label="Substrate S")
plt.plot(sol.t, sol.y[2], label="Product P")
plt.xlabel("Time (h)")
plt.ylabel("Concentration (g/L)")
plt.legend()
plt.show()
