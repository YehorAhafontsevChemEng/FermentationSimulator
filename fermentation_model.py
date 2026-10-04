import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

def get_input(prompt, default_val):
    user_val = input(f"{prompt} [default: {default_val}]: ").strip()
    return float(user_val) if user_val else float(default_val)


print("=== Fermentation Simulator Configuration ===")
print("Press Enter to accept default values for E. coli / Yeast kinetics.\n")

mu_max = get_input("Enter maximum specific growth rate mu_max (1/h)", 0.4)
Ks = get_input("Enter half-saturation constant Ks (g/L)", 0.5)
Yxs = get_input("Enter biomass yield coefficient Yxs (g cells / g sub)", 0.5)
alpha = get_input("Enter growth-associated product yield alpha (g prod / g cells)", 0.3)
beta = get_input("Enter non-growth product rate beta (g prod / g cells / h)", 0.05)
tau = get_input("Enter lag time constant tau in hours (0 = no lag)", 1.0)
Yps = get_input("Enter product yield from substrate Yps (g prod / g sub)", 0.5)

def mu(S):
    return mu_max * S / (Ks + S)

def lag_factor(t):
    if tau == 0:
        return 1.0
    return 1 - np.exp(-t / tau)

def model(t, y):
    X, S, P = y
    S = max(S, 0)
    sat = S / (Ks + S) 

    growth = lag_factor(t) * mu(S)

    dXdt = growth * X
    dPdt = (alpha * growth + beta * sat) * X
    dSdt = -(1 / Yxs) * growth * X - (1 / Yps) * dPdt

    return [dXdt, dSdt, dPdt]


t_span = (0, 24)
t_eval = np.linspace(t_span[0], t_span[1], 300)
y0 = [0.1, 20.0, 0.0]  # [X0, S0, P0] в г/Л

sol = solve_ivp(model, t_span, y0, t_eval=t_eval, rtol=1e-8, atol=1e-10)

plt.figure(figsize=(9, 5))
plt.plot(sol.t, sol.y[0], label="Biomass X (g/L)", linewidth=2)
plt.plot(sol.t, sol.y[1], label="Substrate S (g/L)", linewidth=2, linestyle="--")
plt.plot(sol.t, sol.y[2], label="Product P (g/L)", linewidth=2)

plt.xlabel("Time (h)", fontsize=11)
plt.ylabel("Concentration (g/L)", fontsize=11)
plt.title("Batch Fermentation Simulator (Luedeking-Piret Kinetics)", fontsize=12)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
