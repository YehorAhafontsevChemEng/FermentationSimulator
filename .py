import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

mu_max = float(input("Enter the maximum specific growth rate (mu_max): "))
Ks = float(input("Enter the half-saturation constant (Ks): "))
Yxs = float(input("Enter the yield coefficient (Yxs): "))
alpha = float(input("Enter the product yield coefficient (alpha): "))

def mu(S):
    return mu_max * S / (Ks + S)

def model(t, y):
    X, S, P = y
    S = max(S, 0)           # stops S going slightly negative
    growth = mu(S)
    dXdt = growth * X
    dSdt = -(1 / Yxs) * growth * X
    dPdt = alpha * growth * X
    return [dXdt, dSdt, dPdt]

sol = solve_ivp(model, [0, 24], [0.1, 20, 0],
                t_eval=np.linspace(0, 24, 200))

plt.plot(sol.t, sol.y[0], label="Biomass X")
plt.plot(sol.t, sol.y[1], label="Substrate S")
plt.plot(sol.t, sol.y[2], label="Product P")
plt.xlabel("Time (h)")
plt.ylabel("Concentration (g/L)")
plt.legend()
plt.show()
