"""
experiments.py - Execução automática dos experimentos científicos
Gera dados e figuras para análise e artigo científico.

Metodologia Monte Carlo: cada cenário é executado N_RUNS vezes com sementes
independentes. Os gráficos exibem média ± IC 95% sobre as réplicas.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from simulation import Simulation
from visualization import (plot_epidemic_curves, plot_grid, plot_comparison,
                           plot_comparison_mc, plot_combined)
from config import (DEFAULT_CONFIG, MODEL_SIR, MODEL_SIS, MODEL_SEIR,
                    NEIGHBORHOOD_VON_NEUMANN, NEIGHBORHOOD_MOORE)


OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          '..', 'figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)

N_RUNS = 30  # réplicas Monte Carlo por cenário


# ---------------------------------------------------------------------------
# Utilitário Monte Carlo
# ---------------------------------------------------------------------------

def _mc_sweep(base_cfg: dict, param_key: str, param_values: list,
              metric: str = 'I', label_fmt: str = None) -> dict:
    """
    Varre param_values, executa N_RUNS simulações por valor e retorna
    {label: (steps_arr, mean, ci_lo, ci_hi)} prontos para plot_comparison_mc.
    """
    results_mc = {}
    max_steps = base_cfg.get('steps', 300) + 1

    for val in param_values:
        cfg = {**base_cfg, param_key: val}
        label = label_fmt.format(val) if label_fmt else f'{param_key}={val}'
        runs = []

        for run in range(N_RUNS):
            seed = base_cfg.get('seed', 42) + run * 17
            sim = Simulation({**cfg, 'seed': seed})
            stats = sim.run()
            df = stats.to_dataframe()
            series = df[metric].values
            # Alinha comprimento: preenche com último valor se terminou cedo
            if len(series) < max_steps:
                series = np.pad(series, (0, max_steps - len(series)), mode='edge')
            else:
                series = series[:max_steps]
            runs.append(series)

        arr = np.array(runs)           # (N_RUNS, max_steps)
        mean = arr.mean(axis=0)
        ci = 1.96 * arr.std(axis=0) / np.sqrt(N_RUNS)
        steps_arr = np.arange(max_steps)
        results_mc[label] = (steps_arr, mean, mean - ci, mean + ci)

    return results_mc


def _mc_sweep_summary(base_cfg: dict, param_key: str, param_values: list,
                      label_fmt: str = None) -> pd.DataFrame:
    """
    Gera tabela resumo (pico, passo do pico, total R) com média ± dp
    sobre N_RUNS réplicas.
    """
    rows = []
    for val in param_values:
        cfg = {**base_cfg, param_key: val}
        label = label_fmt.format(val) if label_fmt else f'{param_key}={val}'
        peaks, peak_steps, finals_r = [], [], []

        for run in range(N_RUNS):
            seed = base_cfg.get('seed', 42) + run * 17
            sim = Simulation({**cfg, 'seed': seed})
            stats = sim.run()
            df = stats.to_dataframe()
            peaks.append(df['I'].max())
            peak_steps.append(int(df['I'].idxmax()))
            if 'R' in df.columns:
                finals_r.append(df['R'].iloc[-1])

        def fmt(v):
            return f'{np.mean(v):.1f} ± {np.std(v):.1f}'

        row = {'Parâmetro': label,
               'Pico Infectados': fmt(peaks),
               'Passo do Pico': fmt(peak_steps)}
        if finals_r:
            row['Total Recuperados'] = fmt(finals_r)
        rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Experimentos
# ---------------------------------------------------------------------------

def exp1_beta_variation():
    """Experimento 1: Variação da taxa de transmissão β (γ=0,14 — calibrado influenza)"""
    print(f"  [1/5] Variação de β  ({N_RUNS} réplicas por cenário)...")
    betas = [0.10, 0.18, 0.30, 0.50, 0.70]
    # γ=0,14 ≈ 1/7 dia — duração infecciosa da influenza (Biggerstaff et al., 2014)
    base_cfg = {**DEFAULT_CONFIG, 'gamma': 0.14, 'grid_size': 80, 'steps': 300, 'seed': 42}

    results_mc = _mc_sweep(base_cfg, 'beta', betas, metric='I',
                           label_fmt='β={:.2f}')

    fig = plot_comparison_mc(
        results_mc,
        metric_label='Indivíduos Infectados',
        title=f'Exp. 1 — Influência da Taxa de Transmissão (β)  [n={N_RUNS}]',
        save_path=os.path.join(OUTPUT_DIR, 'exp1_beta_variation.png'))
    plt.close(fig)

    return _mc_sweep_summary(base_cfg, 'beta', betas, label_fmt='β={:.2f}')


def exp2_neighborhood_comparison():
    """Experimento 2: Von Neumann vs Moore"""
    print(f"  [2/5] Comparação de vizinhança  ({N_RUNS} réplicas por cenário)...")
    base_cfg = {**DEFAULT_CONFIG, 'grid_size': 80, 'steps': 300, 'seed': 42}
    neighborhoods = [NEIGHBORHOOD_VON_NEUMANN, NEIGHBORHOOD_MOORE]
    nh_labels = {NEIGHBORHOOD_VON_NEUMANN: 'Von Neumann (4)',
                 NEIGHBORHOOD_MOORE: 'Moore (8)'}

    results_mc = {}
    max_steps = base_cfg['steps'] + 1

    for nh in neighborhoods:
        cfg = {**base_cfg, 'neighborhood': nh}
        label = nh_labels[nh]
        runs = []
        for run in range(N_RUNS):
            seed = base_cfg['seed'] + run * 17
            sim = Simulation({**cfg, 'seed': seed})
            stats = sim.run()
            df = stats.to_dataframe()
            series = df['I'].values
            if len(series) < max_steps:
                series = np.pad(series, (0, max_steps - len(series)), mode='edge')
            else:
                series = series[:max_steps]
            runs.append(series)
        arr = np.array(runs)
        mean = arr.mean(axis=0)
        ci = 1.96 * arr.std(axis=0) / np.sqrt(N_RUNS)
        results_mc[label] = (np.arange(max_steps), mean, mean - ci, mean + ci)

    fig = plot_comparison_mc(
        results_mc,
        title=f'Exp. 2 — Von Neumann vs Moore  [n={N_RUNS}]',
        save_path=os.path.join(OUTPUT_DIR, 'exp2_neighborhood.png'))
    plt.close(fig)

    # Snapshots com semente fixa (apenas para visualização espacial)
    fig2, axes = plt.subplots(2, 4, figsize=(16, 8))
    fig2.suptitle('Exp. 2 — Progressão Espacial: Von Neumann vs Moore', fontsize=13)
    from visualization import CA_CMAP, CA_NORM
    from config import STATE_COLORS, STATE_LABELS
    import matplotlib.patches as mpatches

    for row_idx, nh in enumerate(neighborhoods):
        snap_states, snap_steps_list = [], []

        def collect(step, counts, state, _snap=snap_states, _steps=snap_steps_list):
            if step in [0, 20, 60, 150]:
                _snap.append(state.copy())
                _steps.append(step)

        sim = Simulation({**base_cfg, 'neighborhood': nh})
        sim.run(callback=collect)

        for col_idx, (state, step) in enumerate(zip(snap_states, snap_steps_list)):
            ax = axes[row_idx][col_idx]
            ax.imshow(state, cmap=CA_CMAP, norm=CA_NORM, interpolation='nearest')
            ax.set_title(f'{nh_labels[nh]}\nt={step}', fontsize=8)
            ax.axis('off')

    plt.tight_layout()
    fig2.savefig(os.path.join(OUTPUT_DIR, 'exp2_snapshots.png'), dpi=130, bbox_inches='tight')
    plt.close(fig2)

    rows = []
    for nh in neighborhoods:
        cfg = {**base_cfg, 'neighborhood': nh}
        peaks, durations, finals_r = [], [], []
        for run in range(N_RUNS):
            seed = base_cfg['seed'] + run * 17
            sim = Simulation({**cfg, 'seed': seed})
            stats = sim.run()
            df = stats.to_dataframe()
            peaks.append(df['I'].max())
            dur = df[df['I'] == 0].index[0] if (df['I'] == 0).any() else len(df)
            durations.append(dur)
            finals_r.append(df['R'].iloc[-1])

        def fmt(v):
            return f'{np.mean(v):.1f} ± {np.std(v):.1f}'

        rows.append({'Vizinhança': nh_labels[nh],
                     'Pico I': fmt(peaks),
                     'Duração': fmt(durations),
                     'Total R': fmt(finals_r)})
    return pd.DataFrame(rows)


def exp3_model_comparison():
    """Experimento 3: SIR vs SIS vs SEIR"""
    print(f"  [3/5] Comparação de modelos  ({N_RUNS} réplicas por cenário)...")
    models = [MODEL_SIR, MODEL_SIS, MODEL_SEIR]
    base_cfg = {**DEFAULT_CONFIG, 'grid_size': 80, 'steps': 350, 'seed': 42, 'sigma': 0.08}
    max_steps = base_cfg['steps'] + 1

    results_mc = {}
    for model in models:
        cfg = {**base_cfg, 'model': model}
        runs = []
        for run in range(N_RUNS):
            seed = base_cfg['seed'] + run * 17
            sim = Simulation({**cfg, 'seed': seed})
            stats = sim.run()
            df = stats.to_dataframe()
            series = df['I'].values
            if len(series) < max_steps:
                series = np.pad(series, (0, max_steps - len(series)), mode='edge')
            else:
                series = series[:max_steps]
            runs.append(series)
        arr = np.array(runs)
        mean = arr.mean(axis=0)
        ci = 1.96 * arr.std(axis=0) / np.sqrt(N_RUNS)
        results_mc[model] = (np.arange(max_steps), mean, mean - ci, mean + ci)

    fig_i = plot_comparison_mc(
        results_mc,
        title=f'Exp. 3 — Infectados: SIR vs SIS vs SEIR  [n={N_RUNS}]',
        save_path=os.path.join(OUTPUT_DIR, 'exp3_models_infected.png'))
    plt.close(fig_i)

    # Curvas completas por modelo (semente fixa, apenas ilustração)
    for model in models:
        cfg = {**base_cfg, 'model': model}
        sim = Simulation(cfg)
        stats = sim.run()
        df = stats.to_dataframe()
        fig = plot_epidemic_curves(df, model_type=model,
                                   title=f'Curva Epidemiológica — Modelo {model}')
        fig.savefig(os.path.join(OUTPUT_DIR, f'exp3_curve_{model}.png'),
                    dpi=130, bbox_inches='tight')
        plt.close(fig)

    rows = []
    for model in models:
        cfg = {**base_cfg, 'model': model}
        peaks, finals_s = [], []
        for run in range(N_RUNS):
            seed = base_cfg['seed'] + run * 17
            sim = Simulation({**cfg, 'seed': seed})
            stats = sim.run()
            df = stats.to_dataframe()
            peaks.append(df['I'].max())
            finals_s.append(df['S'].iloc[-1])

        def fmt(v):
            return f'{np.mean(v):.1f} ± {np.std(v):.1f}'

        rows.append({'Modelo': model, 'Pico I': fmt(peaks), 'S Final': fmt(finals_s)})
    return pd.DataFrame(rows)


def exp4_population_density():
    """Experimento 4: Efeito do tamanho da grade (densidade)"""
    print(f"  [4/5] Densidade populacional  ({N_RUNS} réplicas por cenário)...")
    sizes = [40, 60, 80, 100, 120]
    base_cfg = {**DEFAULT_CONFIG, 'initial_infected': 3, 'steps': 300, 'seed': 42}

    results_mc = {}
    for size in sizes:
        cfg = {**base_cfg, 'grid_size': size}
        N = size * size
        label = f'N={size}² ({N})'
        max_steps = cfg['steps'] + 1
        runs = []
        for run in range(N_RUNS):
            seed = base_cfg['seed'] + run * 17
            sim = Simulation({**cfg, 'seed': seed})
            stats = sim.run()
            df = stats.to_dataframe()
            series = df['I'].values / N * 100   # normalizado em %
            if len(series) < max_steps:
                series = np.pad(series, (0, max_steps - len(series)), mode='edge')
            else:
                series = series[:max_steps]
            runs.append(series)
        arr = np.array(runs)
        mean = arr.mean(axis=0)
        ci = 1.96 * arr.std(axis=0) / np.sqrt(N_RUNS)
        results_mc[label] = (np.arange(max_steps), mean, mean - ci, mean + ci)

    fig = plot_comparison_mc(
        results_mc,
        metric_label='Infectados (%)',
        title=f'Exp. 4 — Influência do Tamanho da Grade  [n={N_RUNS}]',
        save_path=os.path.join(OUTPUT_DIR, 'exp4_population.png'))
    plt.close(fig)

    rows = []
    for size in sizes:
        cfg = {**base_cfg, 'grid_size': size}
        N = size * size
        peaks = []
        for run in range(N_RUNS):
            seed = base_cfg['seed'] + run * 17
            sim = Simulation({**cfg, 'seed': seed})
            stats = sim.run()
            df = stats.to_dataframe()
            peaks.append(df['I'].max() / N * 100)
        rows.append({'Grade': f'N={size}² ({N})',
                     'Pico I (%)': f'{np.mean(peaks):.2f} ± {np.std(peaks):.2f}'})
    return pd.DataFrame(rows)


def exp5_initial_infected():
    """Experimento 5: Número inicial de infectados"""
    print(f"  [5/5] Infectados iniciais  ({N_RUNS} réplicas por cenário)...")
    n_inits = [1, 3, 5, 10, 20, 50]
    base_cfg = {**DEFAULT_CONFIG, 'grid_size': 80, 'steps': 300, 'seed': 42}

    results_mc = _mc_sweep(base_cfg, 'initial_infected', n_inits, metric='I',
                           label_fmt='I₀={}')

    fig = plot_comparison_mc(
        results_mc,
        title=f'Exp. 5 — Influência dos Infectados Iniciais (I₀)  [n={N_RUNS}]',
        save_path=os.path.join(OUTPUT_DIR, 'exp5_initial_infected.png'))
    plt.close(fig)

    return _mc_sweep_summary(base_cfg, 'initial_infected', n_inits, label_fmt='I₀={}')


def run_all_experiments():
    """Executa todos os experimentos e salva tabelas."""
    print(f"\n=== Executando Experimentos (Monte Carlo, n={N_RUNS}) ===")
    tables = {}
    tables['beta'] = exp1_beta_variation()
    tables['neighborhood'] = exp2_neighborhood_comparison()
    tables['models'] = exp3_model_comparison()
    tables['population'] = exp4_population_density()
    tables['initial'] = exp5_initial_infected()

    print("\n=== Resumo dos Resultados (média ± dp) ===")
    for name, df in tables.items():
        print(f"\n--- Experimento: {name} ---")
        print(df.to_string(index=False))

    return tables


if __name__ == '__main__':
    run_all_experiments()
    print(f"\nFiguras salvas em: {os.path.abspath(OUTPUT_DIR)}")
