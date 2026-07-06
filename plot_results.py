"""
Plot mean glucose and HbA1c trajectories by diet, for a chosen condition group.

Usage:
    python plot_results.py mean_trajectories_by_condition_120d.csv T2D
    python plot_results.py mean_trajectories_by_condition_365d.csv None

Reads the CSV produced by simulate.py and saves a two-panel PNG:
  top panel    - daily mean glucose, one line per diet
  bottom panel - HbA1c (properly lagged), one line per diet
"""

import sys
import csv
import matplotlib.pyplot as plt

DIETS = ["Western", "Mediterranean", "Plant", "Keto"]
DIET_COLORS = {
    "Western": "#c0392b",
    "Mediterranean": "#2980b9",
    "Plant": "#27ae60",
    "Keto": "#8e44ad",
}


def load_trajectories(path):
    rows = list(csv.DictReader(open(path)))
    days = [int(r['day']) for r in rows]
    return rows, days


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else 'mean_trajectories_by_condition_120d.csv'
    condition = sys.argv[2] if len(sys.argv) > 2 else 'T2D'

    rows, days = load_trajectories(csv_path)

    fig, (ax_g, ax_h) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

    for diet in DIETS:
        g_col = f'{diet}_{condition}_daily_mean_G'
        h_col = f'{diet}_{condition}_H'
        if g_col not in rows[0]:
            print(f'Warning: column {g_col} not found - skipping {diet}')
            continue
        g_vals = [float(r[g_col]) for r in rows]
        h_vals = [float(r[h_col]) for r in rows]
        color = DIET_COLORS.get(diet, None)
        ax_g.plot(days, g_vals, label=diet, color=color, linewidth=2)
        ax_h.plot(days, h_vals, label=diet, color=color, linewidth=2)

    ax_g.set_ylabel('Daily mean glucose (mg/dL)')
    ax_g.set_title(f'Simulated glucose and HbA1c by diet - condition group: {condition}')
    ax_g.legend(loc='best')
    ax_g.grid(alpha=0.3)

    ax_h.set_ylabel('HbA1c (%)')
    ax_h.set_xlabel('Day')
    ax_h.legend(loc='best')
    ax_h.grid(alpha=0.3)

    plt.tight_layout()
    out_name = f'diet_comparison_{condition}.png'
    plt.savefig(out_name, dpi=150)
    print(f'Saved {out_name}')


if __name__ == '__main__':
    main()
