import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# ==============================
# pSEIRS MODEL DEFINITION
# ==============================
def pSEIRS_rhs(y, t, params):
    S, E, I, R = y
    beta, mu, eps, alpha, gamma, omega, tau, p = params
    N = S + E + I + R

    dS = beta * (1 - N) - mu*S - gamma*S*I
    dE = gamma*S*I - mu*E - gamma*S*I*np.exp(-mu*omega)
    dI = gamma*S*I*np.exp(-mu*omega) - (mu + eps + alpha)*I
    dR = p*alpha*I - mu*R
    return np.array([dS, dE, dI, dR])

def rk4_solver(f, y0, t, params):
    n = len(t)
    y = np.zeros((n, len(y0)))
    y[0] = y0
    h = t[1] - t[0]
    for i in range(1, n):
        k1 = f(y[i-1], t[i-1], params)
        k2 = f(y[i-1] + 0.5*h*k1, t[i-1] + 0.5*h, params)
        k3 = f(y[i-1] + 0.5*h*k2, t[i-1] + 0.5*h, params)
        k4 = f(y[i-1] + h*k3, t[i-1] + h, params)
        y[i] = y[i-1] + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    return y

# ==============================
# STREAMLIT UI
# ==============================
st.set_page_config(page_title="pSEIRS Model Simulator", layout="wide")
st.title("🧩 pSEIRS Virus Propagation Simulator")
st.markdown("""
Interactively explore the **probabilistic SEIRS (pSEIRS)** model that describes virus or epidemic spread
with **delays** and **temporary immunity**.
Adjust parameters below to see how infection dynamics change.
""")

# Sidebar parameters
st.sidebar.header("🔧 Model Parameters")
beta = st.sidebar.slider("β (birth rate)", 0.0, 0.5, 0.33, 0.01)
mu = st.sidebar.slider("μ (natural death rate)", 0.0, 0.05, 0.006, 0.001)
eps = st.sidebar.slider("ε (infection death)", 0.0, 0.1, 0.06, 0.005)
alpha = st.sidebar.slider("α (recovery rate)", 0.0, 0.1, 0.04, 0.005)
gamma = st.sidebar.slider("γ (infection rate)", 0.0, 1.0, 0.308, 0.01)
omega = st.sidebar.slider("ω (latency delay)", 0.0, 40.0, 0.15, 0.5)
tau = st.sidebar.slider("τ (immunity delay)", 0.0, 60.0, 30.0, 1.0)
p = st.sidebar.slider("p (immunity probability)", 0.0, 1.0, 1.0, 0.05)
T = st.sidebar.slider("Simulation time (days)", 100, 800, 400, 50)
h = st.sidebar.slider("Step size", 0.001, 0.1, 0.01, 0.001)

# Initial conditions
st.sidebar.header("📈 Initial Fractions")
S0 = st.sidebar.number_input("S(0)", 0.0, 1.0, 0.86)
E0 = st.sidebar.number_input("E(0)", 0.0, 1.0, 0.07)
I0 = st.sidebar.number_input("I(0)", 0.0, 1.0, 0.07)
R0 = st.sidebar.number_input("R(0)", 0.0, 1.0, 0.0)

# Run simulation
params = (beta, mu, eps, alpha, gamma, omega, tau, p)
y0 = np.array([S0, E0, I0, R0])
t = np.arange(0, T, h)
Y = rk4_solver(pSEIRS_rhs, y0, t, params)
S, E, I, R = Y.T

# ==============================
# PLOTS
# ==============================
st.subheader("📊 Evolution of Compartments")
fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(t, S, label="S(t) – Susceptible", color='tab:blue')
ax.plot(t, E, label="E(t) – Exposed", color='tab:orange')
ax.plot(t, I, label="I(t) – Infectious", color='tab:green')
ax.plot(t, R, label="R(t) – Recovered", color='tab:red')
ax.set_xlabel("Time (days)")
ax.set_ylabel("Fraction of population")
ax.set_title("pSEIRS Model Dynamics")
ax.legend()
ax.grid(True)
st.pyplot(fig)

# ==============================
# EXTRA ANALYSIS
# ==============================
st.subheader("📈 Key Observations")
st.markdown(f"""
- **Peak infection (I)**: {I.max():.3f} at day {t[I.argmax()]:.1f}  
- **Final recovered fraction**: {R[-1]:.3f}  
- **Equilibrium population sum**: {S[-1] + E[-1] + I[-1] + R[-1]:.3f}  
""")

st.info("""
🧠 *Tip:* Try reducing `p` to simulate loss of immunity, or increasing `ω` to introduce a longer latency.
Notice how infection persistence and oscillations change dynamically.
""")
