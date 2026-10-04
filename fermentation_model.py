
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares


DEFAULTS = dict(
    mu_max=0.4,   # 1/h    fastest possible growth rate
    Ks=0.05,      # g/L    substrate level where growth is half of max
    Yxs=0.5,      # g/g    g cells made per g substrate
    Yps=0.5,      # g/g    g product made per g substrate
    alpha=0.2,    # g/g    growth-linked product
    beta=0.05,    # g/g/h  non-growth-linked product
    tau=2.0,      # h      lag time constant
    kd=0.04,      # 1/h    death rate when starved
    X0=0.1,       # g/L    starting cells
    S0=4.0,       # g/L    starting substrate
)


def rhs(t, y, p):
    X, S, P = y
    S = max(S, 0.0)
    sat = S / (p["Ks"] + S)                       # ~1 with food, ~0 without
    lag = 1.0 if p["tau"] == 0 else 1 - np.exp(-t / p["tau"])
    growth = lag * p["mu_max"] * sat              # specific growth rate
    death = p["kd"] * (1 - sat)                   # starvation death
    dX = growth * X - death * X
    dP = (p["alpha"] * growth + p["beta"] * sat) * X
    dS = -(1 / p["Yxs"]) * growth * X - (1 / p["Yps"]) * dP
    return [dX, dS, dP]


def simulate(p, t_end=24, n=481):
    t = np.linspace(0, t_end, n)
    sol = solve_ivp(rhs, [0, t_end], [p["X0"], p["S0"], 0.0], args=(p,),
                    t_eval=t, rtol=1e-8, atol=1e-10)
    return sol.t, sol.y[0], sol.y[1], sol.y[2]


# OD600 values read by eye off a published growth curve (approximate!).
# The flat part at 12 h is kept because without it the fit cannot tell how much
# substrate there was (growth that stops needs a plateau to pin down S0). Note
# the original curve is perfectly flat there, which may mean the reading was
# capped, so treat the fitted S0 as approximate.
t_data = np.array([0.0, 3.0, 6.0, 9.0, 12.0])
od_data = np.array([0.3, 0.6, 1.5, 4.0, 4.0])
OD_PER_GL = 1 / 0.35   # OD units per g/L of dry cells (rough rule of thumb)


def fit_to_od():
    """Adjust mu_max, tau and S0 until the simulated OD matches the data."""
    base = dict(DEFAULTS)
    base["X0"] = od_data[0] / OD_PER_GL   # starting cells from the first reading

    def residuals(theta):
        p = dict(base, mu_max=theta[0], tau=theta[1], S0=theta[2])
        t, X, S, P = simulate(p, t_end=t_data[-1], n=361)
        od_model = np.interp(t_data, t, X) * OD_PER_GL
        return np.log(od_model) - np.log(od_data)   # log: equal weight for small and large OD

    res = least_squares(residuals, x0=[0.4, 2.0, 2.8],
                        bounds=([0.05, 0.0, 0.5], [1.5, 6.0, 10.0]))
    fitted = dict(base, mu_max=res.x[0], tau=res.x[1], S0=res.x[2])
    return fitted, res


TURNAROUND = 4.0   # h to empty, clean and refill the vessel between runs


def productivity(t, P):
    """Product per hour of vessel time, counting turnaround between runs."""
    return P / (t + TURNAROUND)


def best_stop_time(p, t_end=36):
    t, X, S, P = simulate(p, t_end=t_end, n=721)
    prod = productivity(t, P)
    i = prod.argmax()
    return t[i], P[i], prod[i], (t, prod)


if __name__ == "__main__":
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))

    # (a) baseline run
    t, X, S, P = simulate(DEFAULTS)
    a = ax[0, 0]
    a.plot(t, X, label="Biomass X")
    a.plot(t, S, label="Substrate S")
    a.plot(t, P, label="Product P")
    a.set(title="A. Baseline batch run", xlabel="Time (h)", ylabel="g/L")
    a.legend()

    # (b) fit to data
    fitted, res = fit_to_od()
    tf, Xf, Sf, Pf = simulate(fitted, t_end=12)
    b = ax[0, 1]
    b.plot(tf, Xf * OD_PER_GL, label="Model (fitted)")
    b.scatter(t_data, od_data, color="black", zorder=3, label="Data (read off figure)")
    b.set(title="B. Fitted to a real OD600 curve", xlabel="Time (h)", ylabel="OD600")
    b.legend()

    # (c) what-if: starting substrate
    c = ax[1, 0]
    for S0 in (2, 4, 8):
        t, X, S, P = simulate(dict(DEFAULTS, S0=S0))
        c.plot(t, P, label=f"S0 = {S0} g/L")
    c.set(title="C. What if we start with more substrate?",
          xlabel="Time (h)", ylabel="Product P (g/L)")
    c.legend()

    # (d) best stop time
    d = ax[1, 1]
    for S0 in (2, 4, 8):
        p = dict(DEFAULTS, S0=S0)
        ts, Pbest, prodbest, (tt, prod) = best_stop_time(p)
        line, = d.plot(tt, prod, label=f"S0 = {S0} g/L")
        d.scatter([ts], [prodbest], color=line.get_color(), zorder=3)
    d.set(title="D. Best time to stop (dots)", xlabel="Stop time (h)",
          ylabel="Product per vessel-hour (g/L/h)")
    d.legend()

    plt.tight_layout()
    plt.savefig("summary_figure.png", dpi=130)

    # printed results
    print("FIT TO DATA")
    print(f"  mu_max = {fitted['mu_max']:.3f} 1/h, lag tau = {fitted['tau']:.2f} h, "
          f"S0 = {fitted['S0']:.2f} g/L")
    print(f"  fit residual (log OD, RMS) = {np.sqrt(np.mean(res.fun**2)):.3f}")
    print("BEST STOP TIME (baseline parameters, 4 h turnaround)")
    for S0 in (2, 4, 8):
        ts, Pb, pr, _ = best_stop_time(dict(DEFAULTS, S0=S0))
        print(f"  S0 = {S0}: stop at {ts:.1f} h, product {Pb:.2f} g/L, "
              f"{pr:.3f} g/L per vessel-hour")
    plt.show()
