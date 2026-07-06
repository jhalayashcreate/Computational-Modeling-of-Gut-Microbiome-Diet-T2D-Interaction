"""
SCFA / T2D multi-layer glucose model - Python implementation.

Reads:
  parameters.csv         - scalar model parameters (name, value)
  diet_targets.csv        - per-diet SCFA/weight steady-state targets
  survey_covariates.csv   - bootstrap pool of real respondent covariates + condition label

Writes:
  mean_trajectories_by_condition.csv  - daily mean of each state variable, by diet x condition
  validation_summary.csv              - simulated H(120) vs condition group, for validation
  sample_patient_intraday.csv         - one example patient's full within-day glucose/insulin
                                         trace on day 60, to sanity-check meal-response shape

Run:  python simulate.py
"""

import csv
import math
import random
from collections import defaultdict, deque

# ---------------------------------------------------------------------------
# LAYER 4: HbA1c convolution kernel
# ---------------------------------------------------------------------------
# Clinical consensus (Labcorp, NGSP, multiple population studies): the most
# recent 30 days contribute ~50% of HbA1c, days 30-60 ~25%, days 60-120 ~25%,
# reflecting the RBC age distribution over the ~120-day RBC lifespan.
# Fitted exponential decay tau=51.2 days reproduces this (~49/27/24).
H_KERNEL_TAU = 51.2
H_KERNEL_WINDOW = 120
H_KERNEL = [math.exp(-a / H_KERNEL_TAU) for a in range(H_KERNEL_WINDOW)]
H_KERNEL_SUM = sum(H_KERNEL)

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
RUNS_PER_DIET = 500
DAYS = 120
DAYS_LONG = 365
DIETS = ["Western", "Mediterranean", "Plant", "Keto"]
CONDITIONS = ["None", "T2D", "Hypertension", "Atherosclerosis", "Other"]

# Meal schedule: (clock_minute, carb_grams) - simple 3-meal day
MEAL_TIMES = [7 * 60, 12 * 60, 18 * 60]  # 7am/12pm/6pm
DT_MIN = None               # set from fast_dt_min parameter at runtime
MINUTES_PER_DAY = 24 * 60

random.seed(42)

# ---------------------------------------------------------------------------
# LOAD INPUTS
# ---------------------------------------------------------------------------
def load_params(path):
    p = {}
    with open(path, newline='') as f:
        for row in csv.DictReader(f):
            p[row['parameter']] = float(row['value'])
    return p


def load_diet_targets(path):
    ss = {}
    with open(path, newline='') as f:
        for row in csv.DictReader(f):
            ss[row['diet']] = {
                'A': float(row['A_ss']), 'P': float(row['P_ss']),
                'B': float(row['B_ss']), 'W': float(row['W_ss']),
            }
    return ss


def load_diet_meals(path):
    meals = {}
    with open(path, newline='') as f:
        for row in csv.DictReader(f):
            meals[row['diet']] = (float(row['carb_per_meal_g']), float(row['fiber_per_meal_g']))
    return meals


def load_covariates(path):
    pool = []
    with open(path, newline='') as f:
        for row in csv.DictReader(f):
            pool.append({
                'bmi': float(row['bmi']),
                'activity': float(row['activity_min_per_week']),
                'stress': float(row['stress_score']),
                'sugar': float(row['sugary_servings']),
                'metformin': row['med_metformin'] == '1',
                'sglt2i': row['med_sglt2i'] == '1',
                'glp1': row['med_glp1'] == '1',
                'insulin': row['med_insulin'] == '1',
                'nonadherence': float(row['med_nonadherence_score']),
                'smoke': row['smoke_current'] == '1',
                'condition': row['condition_label'] if row['condition_label'] in CONDITIONS else 'Other',
            })
    return pool


def clip(x, lo, hi):
    return min(hi, max(lo, x))


# ---------------------------------------------------------------------------
# LAYER 1: intraday Bergman-style meal-response model
# ---------------------------------------------------------------------------
def simulate_intraday(Gb, Ib, M_ratio, beta, p, diet_fiber_factor, carb_per_meal, record=False):
    """
    Simulate one day at fast_dt_min resolution.
    Returns (mean_glucose, peak_glucose, [optional full trace]).
    M_ratio = M(t)/M0_healthy - that day's chronic insulin sensitivity, relative to healthy.
    diet_fiber_factor slows/flattens meal absorption (Layer 2 acute coupling):
      >1 = more fiber = slower, flatter Ra(t).
    """
    p1 = p['p1_SG_healthy'] * (p['p1_SG_floor_ratio'] + (1 - p['p1_SG_floor_ratio']) * M_ratio)
    p2 = p['p2_Bergman']
    p3 = p['p3_base'] * M_ratio
    n = p['n_clearance']
    gamma = p['gamma_pancreas']
    h = Gb  # secretion threshold ~ basal glucose
    k_abs = p['k_abs_base'] / diet_fiber_factor  # slower absorption -> lower k_abs
    Vg = p['Vg']  # whole-body glucose distribution volume, dL

    G, X, I = Gb, 0.0, Ib
    trace = []
    glucose_samples = []

    n_steps = MINUTES_PER_DAY // DT_MIN
    for step in range(n_steps):
        t_min = step * DT_MIN

        # sum glucose appearance from all meals still active
        Ra = 0.0
        for meal_t in MEAL_TIMES:
            dt_since = t_min - meal_t
            if dt_since >= 0:
                Ra += (k_abs ** 2) * dt_since * math.exp(-k_abs * dt_since) * carb_per_meal * p['f_carb_bioavail']
                # f_carb_bioavail: mg systemic glucose per gram ingested carb (calibrated, not 1:1 -
                # accounts for first-pass hepatic uptake and non-glucose digestion products)

        dG = (-p1 * (G - Gb) - X * G + Ra / Vg) * DT_MIN
        dX = (-p2 * X + p3 * (I - Ib)) * DT_MIN
        secretion = gamma * beta * max(G - h, 0.0)
        dI = (secretion - n * (I - Ib)) * DT_MIN

        G = max(G + dG, 20.0)   # floor to avoid runaway negative glucose
        X = X + dX
        I = max(I + dI, 0.0)

        glucose_samples.append(G)
        if record:
            trace.append((t_min, G, I, X))

    mean_g = sum(glucose_samples) / len(glucose_samples)
    peak_g = max(glucose_samples)
    return mean_g, peak_g, trace


# ---------------------------------------------------------------------------
# MAIN SIMULATION
# ---------------------------------------------------------------------------
def run_simulation(days, label):
    p = load_params('parameters.csv')
    global DT_MIN
    DT_MIN = int(p['fast_dt_min'])
    ss = load_diet_targets('diet_targets.csv')
    diet_meals = load_diet_meals('diet_meals.csv')
    pool = load_covariates('survey_covariates.csv')

    S_base = p['beta_A'] + p['beta_P'] + p['beta_B']

    series_names = ["A", "P", "B", "W", "S", "M0", "M", "beta0", "beta", "daily_mean_G", "daily_peak_G", "H_ss", "H"]
    acc = {d: {c: {'sums': {s: [0.0] * (days + 1) for s in series_names}, 'n': 0} for c in CONDITIONS}
           for d in DIETS}

    sample_trace_rows = []
    sample_taken = False

    for d in DIETS:
        targets = ss[d]
        carb_per_meal, fiber_per_meal = diet_meals[d]
        daily_fiber = 3 * fiber_per_meal  # 3 meals/day
        for run in range(RUNS_PER_DIET):
            cov = random.choice(pool)
            cond = cov['condition']

            fiber_factor = 0.8 + 0.4 * random.random()  # per-patient variability in fermentation efficiency
            w_factor = 0.7 + 0.6 * random.random()
            patient_daily_fiber = daily_fiber * fiber_factor
            W_target = targets['W'] * w_factor

            on_any_diabetes_med = cov['metformin'] or cov['sglt2i'] or cov['glp1'] or cov['insulin']
            M0_patient = clip(
                p['M0_healthy'] - p['k_BMI'] * (cov['bmi'] - p['BMI_ref'])
                - (p['k_med'] if on_any_diabetes_med else 0),
                0.5, 10
            )
            G_chronic = clip(
                p['G_healthy_ref'] - p['gamma_chronic'] * (M0_patient - p['M0_healthy']),
                p['G_min'], p['G_max']
            )
            Ib_patient = p['Ib_basal']
            beta0_patient = clip(
                1.0 - p['k_beta0'] * max(G_chronic - p['G_beta_toxic_thresh'], 0),
                p['beta_min'], 1.0
            )

            A, P_, B, W = 1.0, 1.0, 1.0, 0.0
            M = M0_patient

            # Warm-up: run one throwaway intraday day at this patient's starting state to get
            # a real (meal-inclusive) daily mean glucose, rather than pre-filling history with
            # the fasting-level G_chronic - that mismatch was creating an artificial ~120-day
            # ramp-up in every trajectory (H_ss vs H didn't converge until day 120 instead of day 0).
            warmup_mean_g, _, _ = simulate_intraday(
                G_chronic, Ib_patient, M0_patient / p['M0_healthy'], beta0_patient, p,
                diet_fiber_factor=1.0, carb_per_meal=carb_per_meal, record=False
            )
            glucose_history = deque([warmup_mean_g] * H_KERNEL_WINDOW, maxlen=H_KERNEL_WINDOW)
            beta_cell = beta0_patient

            bucket = acc[d][cond]
            bucket['n'] += 1

            for t in range(days + 1):
                S = p['beta_A'] * A + p['beta_P'] * P_ + p['beta_B'] * B + p['beta_W'] * W

                pAdherent = 1 - cov['nonadherence'] / 4
                adherent_today = random.random() < clip(pAdherent, 0, 1)

                bracket = (1
                    + p['beta_S'] * (S / S_base - 1)
                    + p['beta_E'] * (cov['activity'] - p['E_ref']) / p['E_ref']
                    - p['beta_Str'] * (cov['stress'] - 2) / 4
                    - (p['delta_smoke'] if cov['smoke'] else 0)
                    - p['beta_Sugar'] * (cov['sugar'] - p['Sugar_ref']) / (p['Sugar_ref'] + 1)
                    + (p['kappa_metformin'] if (adherent_today and cov['metformin']) else 0)
                    + (p['kappa_sglt2i'] if (adherent_today and cov['sglt2i']) else 0)
                    + (p['kappa_glp1'] if (adherent_today and cov['glp1']) else 0)
                )
                M_target = clip(M0_patient * bracket, 0.3, 12)
                M += (M_target - M) / p['tau_M']

                # ---- LAYER 1: intraday meal-response simulation ----
                M_ratio = M / p['M0_healthy']
                diet_fiber_factor = clip(S / S_base, 0.5, 2.0)  # more SCFA signal -> slower absorption
                record_today = (not sample_taken) and (d == "Western") and (cond == "T2D") and (t == 60)
                mean_g, peak_g, trace = simulate_intraday(
                    G_chronic, Ib_patient, M_ratio, beta_cell, p, diet_fiber_factor,
                    carb_per_meal=carb_per_meal, record=record_today
                )
                if record_today:
                    for (t_min, G_val, I_val, X_val) in trace:
                        sample_trace_rows.append([t_min, G_val, I_val, X_val])
                    sample_taken = True

                if cov['insulin'] and adherent_today:
                    mean_g = max(mean_g - p['kappa_insulin'], p['G_min'])
                    peak_g = max(peak_g - p['kappa_insulin'], p['G_min'])

                Hss = (mean_g + 46.7) / 28.7  # instantaneous ADAG-equivalent, kept for comparison only

                # ---- LAYER 4: HbA1c as a weighted convolution over the real RBC-lifespan window ----
                glucose_history.append(mean_g)
                weighted_sum = 0.0
                for age, g in enumerate(reversed(glucose_history)):
                    weighted_sum += H_KERNEL[age] * g
                mean_weighted_glucose = weighted_sum / H_KERNEL_SUM
                H = (mean_weighted_glucose + 46.7) / 28.7

                # ---- LAYER 3: beta-cell mass update, driven by TODAY's real simulated glucose ----
                regen = p['r_beta_regen'] * max(p['G_beta_low_thresh'] - mean_g, 0) / p['G_beta_low_thresh']
                damage = p['r_beta_damage'] * max(mean_g - p['G_beta_toxic_thresh'], 0) ** 2
                beta_cell = clip(beta_cell * (1 + regen - damage), p['beta_min'], 1.0)

                bucket['sums']['A'][t] += A
                bucket['sums']['P'][t] += P_
                bucket['sums']['B'][t] += B
                bucket['sums']['W'][t] += W
                bucket['sums']['S'][t] += S
                bucket['sums']['M0'][t] += M0_patient
                bucket['sums']['M'][t] += M
                bucket['sums']['beta0'][t] += beta0_patient
                bucket['sums']['beta'][t] += beta_cell
                bucket['sums']['daily_mean_G'][t] += mean_g
                bucket['sums']['daily_peak_G'][t] += peak_g
                bucket['sums']['H_ss'][t] += Hss
                bucket['sums']['H'][t] += H

                A += p['k_A_yield'] * patient_daily_fiber - A / p['tau_A']
                P_ += p['k_P_yield'] * patient_daily_fiber - P_ / p['tau_P']
                B += p['k_B_yield'] * patient_daily_fiber - B / p['tau_B']
                W += (W_target - W) / p['tau_W']

            if (run + 1) % 100 == 0:
                print(f'{d}: {run+1}/{RUNS_PER_DIET} virtual patients done')

    # ---- WRITE mean_trajectories_by_condition CSV ----
    traj_path = f'mean_trajectories_by_condition_{label}.csv'
    with open(traj_path, 'w', newline='') as f:
        w = csv.writer(f)
        header = ['day']
        for d in DIETS:
            for c in CONDITIONS:
                for s in series_names:
                    header.append(f'{d}_{c}_{s}')
        w.writerow(header)
        for t in range(days + 1):
            row = [t]
            for d in DIETS:
                for c in CONDITIONS:
                    n = acc[d][c]['n'] or 1
                    for s in series_names:
                        row.append(acc[d][c]['sums'][s][t] / n)
            w.writerow(row)

    # ---- WRITE validation_summary CSV ----
    val_path = f'validation_summary_{label}.csv'
    with open(val_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['diet', 'condition', f'simulated_H{days}_mean', 'n_virtual_patients'])
        for d in DIETS:
            for c in CONDITIONS:
                n = acc[d][c]['n']
                h_final = acc[d][c]['sums']['H'][days] / n if n else float('nan')
                w.writerow([d, c, h_final, n])

    # ---- WRITE sample_patient_intraday CSV (only meaningful once, but harmless to rewrite) ----
    sample_path = f'sample_patient_intraday_{label}.csv'
    with open(sample_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['minute_of_day', 'glucose_mgdl', 'insulin_uUmL', 'remote_insulin_action'])
        for row in sample_trace_rows:
            w.writerow(row)

    print(f'\n[{label}] Done. Wrote {traj_path}, {val_path}, {sample_path}')


if __name__ == '__main__':
    run_simulation(DAYS, f'{DAYS}d')
    run_simulation(DAYS_LONG, f'{DAYS_LONG}d')
