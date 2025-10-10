import matplotlib.pyplot as plt
import re
import math
import os

def load_points(filename):
    points = []
    with open(filename, 'r') as f:
        for line in f:
            match = re.match(r"x:\s*([-0-9.]+),\s*y:\s*([-0-9.]+),\s*inside_circle:\s*(yes|no)", line.strip())
            if match:
                x = float(match.group(1))
                y = float(match.group(2))
                inside = match.group(3) == "yes"
                points.append((x, y, inside))
    return points

def plot_points(points, title, save_png=False):
    x_in = [p[0] for p in points if p[2]]
    y_in = [p[1] for p in points if p[2]]
    x_out = [p[0] for p in points if not p[2]]
    y_out = [p[1] for p in points if not p[2]]

    fig, ax = plt.subplots(figsize=(6,6))
    ax.scatter(x_in, y_in, s=8, label='Wewnątrz koła')
    ax.scatter(x_out, y_out, s=8, label='Poza kołem')

    circle = plt.Circle((0,0), 1.0, color='blue', fill=False, linewidth=2, label='Koło (r=1)')
    ax.add_artist(circle)

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal', 'box')
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend(loc='upper right')
    plt.grid(True, linestyle='--', alpha=0.5)

    if save_png:
        fname = f"points_plot_{title.replace(' ','_')}.png"
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        print(f"Zapisano: {fname}")

    plt.show()

def _float_after_colon(line):
    """Zwraca float występujący po dwukropku w linii (jeśli jest)."""
    m = re.search(r":\s*([-+]?\d*\.?\d+)", line)
    if m:
        return float(m.group(1))
    return None

def load_results(filename):
    results = {}
    with open(filename, 'r') as f:
        for line in f:
            if "Number of points" in line:
                m = re.search(r":\s*(\d+)", line)
                if m: results["N"] = int(m.group(1))
            elif "Estimated area" in line:
                val = _float_after_colon(line)
                if val is not None: results["estimated"] = val
            elif "Exact area" in line:
                val = _float_after_colon(line)
                if val is not None: results["exact"] = val
            elif "Relative error" in line:
                val = _float_after_colon(line)
                if val is not None: results["error"] = val
    return results

def plot_area_comparison(results_list, save_png=False):
    Ns = [r["N"] for r in results_list]
    estimated = [r["estimated"] for r in results_list]
    exact = [r["exact"] for r in results_list]
    errors = [r.get("error", float('nan')) for r in results_list]

    fig, ax = plt.subplots(figsize=(8,5))
    bar_width = 0.35
    indices = list(range(len(Ns)))

    ax.bar([i - bar_width/2 for i in indices], estimated, width=bar_width, label='Oszacowane pole')
    ax.bar([i + bar_width/2 for i in indices], exact, width=bar_width, label='Rzeczywiste pole πr²')

    for i in indices:
        top = max(estimated[i], exact[i])
        ax.text(i, top + 0.08, f"Błąd: {errors[i]:.2f}%", ha='center', fontsize=9)

    ax.set_xticks(indices)
    ax.set_xticklabels([str(N) for N in Ns])
    ax.set_ylabel("Pole koła")
    ax.set_xlabel("Liczba punktów N")
    ax.set_title("Porównanie oszacowanego i rzeczywistego pola koła")
    ax.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.6)

    if save_png:
        plt.savefig("area_comparison.png", dpi=150, bbox_inches='tight')
        print("Zapisano: area_comparison.png")

    plt.show()

if __name__ == "__main__":
    # Pliki generowane przez C++
    point_files = ["points_1000.txt", "points_10000.txt", "points_100000.txt"]
    for fname in point_files:
        if os.path.exists(fname):
            pts = load_points(fname)
            plot_points(pts, f"Punkty Monte Carlo – {fname}", save_png=False)
        else:
            print(f"Brak pliku: {fname} (pomiń).")

    result_files = ["results_1000.txt", "results_10000.txt", "results_100000.txt"]
    results_data = []
    for fname in result_files:
        if os.path.exists(fname):
            results_data.append(load_results(fname))
        else:
            print(f"Brak pliku: {fname} (pomiń).")

    if results_data:
        plot_area_comparison(results_data, save_png=False)
    else:
        print("Brak plików results_*.txt — uruchom program C++ najpierw.")
