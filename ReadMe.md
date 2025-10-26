# 🧩 pSEIRS Virus Propagation Simulator

Interactive **Streamlit web app** to simulate the *probabilistic SEIRS (pSEIRS)* model of virus or epidemic propagation with **latency** and **temporary immunity**.

---

## 🚀 Features

- Real-time simulation of **S–E–I–R** compartments using a Runge–Kutta 4th order (RK4) integrator.
- Adjustable parameters:
  - Infection, recovery, death, and immunity rates.
  - Probabilistic immunity (`p`) and latency (`ω`) delays.
  - Initial compartment fractions and simulation duration.
- Instant plot of compartment evolution.
- Dynamic key metrics: infection peak, recovery fraction, and equilibrium sum.

---

## 🧠 Model Overview

The **pSEIRS model** extends the classical SEIR epidemiological model by including:
- **Delayed transitions** (latency `ω`, immunity `τ`),
- **Temporary immunity** (parameter `p`),
- And **probabilistic reinfection**.

It’s inspired by:
> M.S.S. Khan (2014). *A Computer Virus Propagation Model Using Delay Differential Equations with Probabilistic Contagion and Immunity.*

---

## ⚙️ How to run locally

1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/pSEIRS-simulator.git
   cd pSEIRS-simulator
