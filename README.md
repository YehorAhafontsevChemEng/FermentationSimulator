# Fermentation Simulator

A Python model of batch fermentation. It simulates cell growth, substrate use and product formation, fits some of the parameters to a real OD600 growth curve, and uses the model to find a sensible harvest time.

![Summary figure](summary_figure.png)

## What it does

1. **Simulates a batch culture** with three variables: biomass (X), substrate (S) and product (P) 
2. Units: all in g/L.
3. **Fits three parameters** (`mu_max`, lag `tau`, starting substrate `S0`) to an OD600 growth curve using `scipy.optimize.least_squares` with logarithmic residuals.
4. **Runs a what-if comparison** of different starting substrate amounts.
5. **Finds the best harvest time** by maximising product per vessel-hour, including a fixed turnaround time between runs.

## Model

Kinetic ingredients:

* **Monod growth:** growth slows as substrate runs low (half-saturation constant $K_s$).
* **Lag phase:** growth ramps up smoothly after inoculation (time constant $\tau$).
* **Luedeking-Piret product formation:** one part tied to growth ($\alpha$), one part not ($\beta$).
* **Starvation death:** cells die at rate $k_d$ once substrate is gone.

The system is three coupled ODEs, solved with `scipy.integrate.solve_ivp`:

$$
\begin{aligned}
\text{Sat} &= \frac{S}{K_s + S} \\
\text{Lag} &= 1 - e^{-t/\tau} \\
\mu &= \text{Lag} \cdot \mu_{max} \cdot \text{Sat} \\
\text{Death} &= k_d \cdot (1 - \text{Sat})
\end{aligned}
$$

1. **Biomass:**
   $$\frac{dX}{dt} = (\mu - \text{Death}) \cdot X$$

2. **Product:**
   $$\frac{dP}{dt} = (\alpha \cdot \mu + \beta \cdot \text{Sat}) \cdot X$$

3. **Substrate** (used to make new cells and product; dying cells do not return substrate):
   $$\frac{dS}{dt} = -\frac{1}{Y_{x/s}} \cdot \mu \cdot X - \frac{1}{Y_{p/s}} \cdot \frac{dP}{dt}$$

## Results

Fitted to the OD600 data below:

| Parameter | Fitted value |
|---|---|
| Max growth rate, `mu_max` | about 0.33 per hour |
| Lag time constant, `tau` | about 0.9 h |
| Starting substrate, `S0` | about 4 g/L |

Harvest-time analysis (baseline parameters, 4 h turnaround): the best moment to stop is when the substrate runs out, around 9 h for `S0` = 4 g/L. After that, product has levelled off and extra time only reduces product per vessel-hour.

## Data

OD600 values were read by eye off a published growth curve of a recombinant *E. coli* batch culture, so they are approximate.

Source: *[add the paper or page the figure came from]*

| Time (h) | OD600 |
|---|---|
| 0 | 0.3 |
| 3 | 0.6 |
| 6 | 1.5 |
| 9 | 4.0 |
| 12 | 4.0 |

## Assumptions and limitations

* **Constant conditions:** pH, temperature and oxygen are assumed constant and are not modelled.
* **Simple kinetics:** growth follows Monod kinetics, and product follows Luedeking-Piret. Product toxicity and byproducts are not included.
* **Batch only:** nothing is added or removed during the run (no fed-batch or continuous operation).
* **Partial fitting:** only `mu_max`, `tau` and `S0` were fitted. `Ks`, the yields (`Yxs`, `Yps`), `alpha`, `beta` and `kd` are placeholder values, so product predictions are illustrative.
* **Calibration, not validation:** the model was fitted to one dataset and has not been tested on a different one.
* **Approximate data:** the points were read off a figure by eye. The flat plateau at 12 h may be a capped reading.
* **OD conversion:** OD is converted to g/L with a rough factor (about 0.35 g/L per OD unit), which depends on strain and instrument.
* **Simple economics:** the harvest-time analysis maximises product per vessel-hour with a fixed turnaround time. It does not include costs, prices or downstream processing.

## Getting started

### Prerequisites

* Python 3.8+
* `numpy`, `scipy`, `matplotlib`

### Installation

```bash
git clone https://github.com/YehorAhafontsevChemEng/FermentationSimulator.git
cd FermentationSimulator
pip install numpy scipy matplotlib
```

### Usage

```bash
python fermentation_sim.py
```

This prints the fitted parameters and best stop times, and saves `summary_figure.png`.

## Possible next steps

* Fed-batch operation (substrate feed and changing volume)
* Product inhibition of growth
* Validating the fit on a second, independent dataset
