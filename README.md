# FermentationSimulator
# Bioprocess Kinetic Simulator & Process Optimizer

A Python-based simulation engine for batch fermentation kinetics, parameter identification, and economic yield optimization. Designed to bridge bioprocess engineering principles with numerical simulation (`SciPy`).

![Simulation Summary](summary_figure.png)

## Key Features

* **Advanced Kinetic Modeling:** 
  * **Monod Substrate-Limited Growth:** Accounts for substrate saturation effects ($K_s$).
  * **Luedeking-Piret Product Formation:** Separates growth-associated ($\alpha$) and non-growth-associated ($\beta$) synthesis.
  * **Lag-Phase Dynamics:** Smooth exponential transition ($\tau$) upon inoculation.
  * **Starvation & Cell Death:** Models biomass decay ($k_d$) upon substrate depletion.
* **Parameter Estimation:** Fits non-linear kinetic parameters ($\mu_{max}$, $\tau$, $S_0$) to experimental $\text{OD}_{600}$ growth data using logarithmic residuals via `scipy.optimize.least_squares`.
* **Process Economics & Optimization:** Evaluates vessel productivity ($\text{g/L/h}$) taking equipment turnaround time (cleaning, sterilization) into account to determine the optimal harvest point.
* ## Experimental Data Source

The parameter estimation model is calibrated against published experimental growth data for recombinant *Escherichia coli*:
* **Source:** *Batch fermentation growth curve of recombinant E. coli clone* ($\text{OD}_{600}$ vs Time).
* **Validation:** Used to estimate maximum specific growth rate ($\mu_{max}$), lag phase constant ($\tau$), and initial substrate concentration ($S_0$) via non-linear least squares optimization.

---

## Mathematical Architecture

The system is governed by three coupled Non-Linear Ordinary Differential Equations (ODEs):

$$\begin{aligned} \text{Sat} &= \frac{S}{K_s + S} \\ \text{Lag} &= 1 - e^{-t/\tau} \\ \mu &= \text{Lag} \cdot \mu_{max} \cdot \text{Sat} \\ \text{Death} &= k_d \cdot (1 - \text{Sat}) \end{aligned}$$

1. **Biomass Accumulation ($X$):**
   $$\frac{dX}{dt} = (\mu - \text{Death}) \cdot X$$

2. **Product Synthesis ($P$):**
   $$\frac{dP}{dt} = (\alpha \cdot \mu + \beta \cdot \text{Sat}) \cdot X$$

3. **Substrate Consumption ($S$):**
   $$\frac{dS}{dt} = -\frac{1}{Y_{x/s}} \left(\frac{dX}{dt}\right) - \frac{1}{Y_{p/s}} \left(\frac{dP}{dt}\right)$$

---

## Getting Started

### Prerequisites

* Python 3.8+
* `numpy`
* `scipy`
* `matplotlib`

### Installation

```bash
git clone [https://github.com/YehorAhafontsevChemEng/FermentationSimulator.git](https://github.com/YehorAhafontsevChemEng/FermentationSimulator.git)
cd FermentationSimulator
pip install numpy scipy matplotlib

## Assumptions and limitations

- **Constant conditions:** pH, temperature and oxygen are assumed constant and are not modelled.
- **Simple kinetics:** growth follows Monod kinetics, and product follows Luedeking-Piret. Product toxicity and byproducts are not included.
- **Batch only:** nothing is added or removed during the run (no fed-batch or continuous operation).
- **Partial fitting:** only `mu_max`, `lag tau` and `S0` were fitted to data. `Ks`, the yields (`Yxs`, `Yps`), `alpha`, `beta` and `kd` are placeholder values, so product predictions are illustrative.
- **Approximate data:** the OD600 points were read by eye off a published figure. The flat plateau at 12 h may be a capped reading.
- **OD conversion:** OD is converted to g/L with a rough factor (about 0.35 g/L per OD unit), which depends on strain and instrument.
- **Simple economics:** the "best stop time" maximises product per vessel-hour, including a fixed turnaround time. It does not include costs, prices or downstream processing.
