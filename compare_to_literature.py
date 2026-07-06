"""
Compare simulated HbA1c change per diet against published RCT/meta-analysis
effect sizes, using a one-sample z-test.

For each published study: SE is back-calculated from its reported 95% CI
    SE = (CI_upper - CI_low) / 3.92
Then: z = (model_delta_H - published_delta_H) / SE
      p = two-tailed p-value for that z

A non-significant p (e.g. p > 0.05) means your model's predicted change is
statistically consistent with the published result - not proof the model is
"correct", but evidence it isn't producing an implausible effect size.

Usage:
    python compare_to_literature.py mean_trajectories_by_condition_120d.csv T2D
    python compare_to_literature.py mean_trajectories_by_condition_365d.csv T2D

Note on duration: published studies vary in length (12 weeks to 12 months).
This script compares against whatever final day is in your CSV (120 or 365)
and prints each study's approximate duration alongside so you can judge
whether the comparison is apples-to-apples or a rougher directional check.
"""

import sys
import csv
import math

DIETS = ["Western", "Mediterranean", "Plant", "Keto"]

# Published effect sizes: HbA1c change (percentage points) vs. a control diet,
# in adults with T2D. ci_low/ci_high are the reported 95% CI bounds.
# duration_days is approximate, drawn from each study's stated trial length.
STUDIES = {
    "Mediterranean": [
        {"source": "Mediterranean meta-analysis, 11 RCTs (PMC12735701, 2025)",
         "delta_h": -0.307, "ci_low": -0.451, "ci_high": -0.163, "duration_days": 168},
        {"source": "Mediterranean meta-analysis, 7 RCTs, n=1371 (BMC Nutrition, 2024)",
         "delta_h": -0.39, "ci_low": -0.58, "ci_high": -0.20, "duration_days": 180},
    ],
    "Plant": [
        {"source": "Vegan diet meta-analysis, 11 trials, n=796, >=12wk (PMC9540559)",
         "delta_h": -0.18, "ci_low": -0.29, "ci_high": -0.07, "duration_days": 84},
        {"source": "Vegetarian/vegan meta-analysis, 8 trials, n=369",
         "delta_h": -0.29, "ci_low": -0.45, "ci_high": -0.12, "duration_days": 126},
        {"source": "Vegetarian/vegan meta-analysis (PubMed 40037300)",
         "delta_h": -0.36, "ci_low": -0.54, "ci_high": -0.19, "duration_days": 126},
    ],
    "Keto": [
        {"source": "Network meta-analysis, ketogenic vs control (PMC10384204)",
         "delta_h": -0.73, "ci_low": -1.19, "ci_high": -0.28, "duration_days": 180},
        {"source": "Network meta-analysis, low-carb vs control (PMC10384204)",
         "delta_h": -0.69, "ci_low": -1.32, "ci_high": -0.06, "duration_days": 180},
        {"source": "Ketogenic diet meta-analysis (PMC9246466) - no CI reported, shown for reference only",
         "delta_h": -1.45, "ci_low": None, "ci_high": None, "duration_days": 180},
    ],
    "Western": [
        {"source": "Conventional diabetic diet control arm, Lee et al. 2016 brown-rice-vegan RCT "
                   "(PLOS One) - NOTE: 'conventional diet' = Korean Diabetes Association standard "
                   "diabetic diet, not an unmodified Western diet; shown as the closest available "
                   "standard-care comparator, not a precise match",
         "delta_h": -0.20, "ci_low": None, "ci_high": None, "duration_days": 84},
    ],
}


def normal_cdf(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def z_test(model_delta, published_delta, ci_low, ci_high):
    se = (ci_high - ci_low) / 3.92
    z = (model_delta - published_delta) / se
    p = 2 * (1 - normal_cdf(abs(z)))
    return z, p, se


def load_delta_h(csv_path, condition):
    rows = list(csv.DictReader(open(csv_path)))
    day0, day_final = rows[0], rows[-1]
    final_day = int(day_final['day'])
    deltas = {}
    for diet in DIETS:
        col = f'{diet}_{condition}_H'
        if col not in day0:
            continue
        h0 = float(day0[col])
        hf = float(day_final[col])
        deltas[diet] = hf - h0
    return deltas, final_day


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else 'mean_trajectories_by_condition_120d.csv'
    condition = sys.argv[2] if len(sys.argv) > 2 else 'T2D'

    deltas, final_day = load_delta_h(csv_path, condition)

    print(f'\nModel run: {csv_path}  |  condition group: {condition}  |  horizon: day 0-{final_day}\n')
    print(f'{"Diet":14s} {"Model dH":>9s}  {"Study dH":>9s}  {"95% CI":>18s}  {"z":>6s}  {"p":>7s}  {"~duration":>10s}  Source')
    print('-' * 130)

    for diet in DIETS:
        model_delta = deltas.get(diet)
        if model_delta is None:
            continue
        studies = STUDIES.get(diet, [])
        if not studies:
            print(f'{diet:14s} {model_delta:9.3f}   (no published comparator - this diet is typically the control arm)')
            continue
        for s in studies:
            if s['ci_low'] is None:
                print(f'{diet:14s} {model_delta:9.3f}  {s["delta_h"]:9.3f}  {"n/a (no CI)":>18s}  {"-":>6s}  {"-":>7s}  '
                      f'{s["duration_days"]:8d}d  {s["source"]}')
                continue
            z, p, se = z_test(model_delta, s['delta_h'], s['ci_low'], s['ci_high'])
            ci_str = f'[{s["ci_low"]:.3f}, {s["ci_high"]:.3f}]'
            flag = '' if p > 0.05 else '  <-- differs significantly'
            print(f'{diet:14s} {model_delta:9.3f}  {s["delta_h"]:9.3f}  {ci_str:>18s}  {z:6.2f}  {p:7.4f}  '
                  f'{s["duration_days"]:8d}d  {s["source"]}{flag}')
        print()

    print('Note: duration_days is each study\'s approximate reported trial length, shown so you can')
    print('judge whether comparing against your day-{} run is a fair like-for-like comparison,'.format(final_day))
    print('or more of a directional sanity check when durations differ substantially.')


if __name__ == '__main__':
    main()
