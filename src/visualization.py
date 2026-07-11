"""
visualization.py - Visualização do autômato celular e curvas epidemiológicas
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.gridspec import GridSpec
import os

from config import S, E, I, R, STATE_COLORS, STATE_LABELS


# Colormap para a grade
_COLORS = [STATE_COLORS[S], STATE_COLORS[E], STATE_COLORS[I], STATE_COLORS[R]]
CA_CMAP = ListedColormap(_COLORS)
CA_NORM = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], CA_CMAP.N)


def plot_grid(state: np.ndarray, step: int = 0, ax=None, title: str = None):
    """Plota o mapa bidimensional do autômato celular."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = ax.get_figure()

    ax.imshow(state, cmap=CA_CMAP, norm=CA_NORM, interpolation='nearest')
    ax.set_title(title or f'Grade CA — Passo {step}', fontsize=12)
    ax.axis('off')

    # Legenda
    patches = [
        mpatches.Patch(color=STATE_COLORS[S], label=STATE_LABELS[S]),
        mpatches.Patch(color=STATE_COLORS[E], label=STATE_LABELS[E]),
        mpatches.Patch(color=STATE_COLORS[I], label=STATE_LABELS[I]),
        mpatches.Patch(color=STATE_COLORS[R], label=STATE_LABELS[R]),
    ]
    ax.legend(handles=patches, loc='lower right', fontsize=8,
              framealpha=0.8, edgecolor='gray')
    return fig


def plot_epidemic_curves(df, model_type: str = 'SIR', ax=None, title: str = None):
    """Plota curvas epidemiológicas S, E, I, R ao longo do tempo."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = ax.get_figure()

    total = df['S'].iloc[0] + df['I'].iloc[0] + df.get('E', 0).iloc[0] + df.get('R', 0).iloc[0]

    ax.plot(df['step'], df['S'] / total * 100, color=STATE_COLORS[S],
            linewidth=2, label='Suscetível (S)')

    if model_type == 'SEIR' and 'E' in df.columns and df['E'].max() > 0:
        ax.plot(df['step'], df['E'] / total * 100, color=STATE_COLORS[E],
                linewidth=2, label='Exposto (E)', linestyle='--')

    ax.plot(df['step'], df['I'] / total * 100, color=STATE_COLORS[I],
            linewidth=2.5, label='Infectado (I)')

    if model_type in ('SIR', 'SEIR'):
        ax.plot(df['step'], df['R'] / total * 100, color=STATE_COLORS[R],
                linewidth=2, label='Recuperado (R)')

    ax.set_xlabel('Passo de Simulação (t)', fontsize=11)
    ax.set_ylabel('Proporção da População (%)', fontsize=11)
    ax.set_title(title or f'Curva Epidemiológica — Modelo {model_type}', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 105)
    return fig


def plot_combined(state: np.ndarray, df, step: int, model_type: str,
                  config: dict, save_path: str = None):
    """Cria figura combinada: grade CA + curvas epidemiológicas."""
    fig = plt.figure(figsize=(14, 6))
    fig.patch.set_facecolor('#f8f8f8')
    gs = GridSpec(1, 2, figure=fig, wspace=0.3)

    ax1 = fig.add_subplot(gs[0, 0])
    plot_grid(state, step=step, ax=ax1)

    ax2 = fig.add_subplot(gs[0, 1])
    plot_epidemic_curves(df, model_type=model_type, ax=ax2)

    # Info de parâmetros
    info = (f"β={config.get('beta', '?')}  γ={config.get('gamma', '?')}  "
            f"N={config.get('grid_size', '?')}²  "
            f"Viz={config.get('neighborhood', '?')}")
    fig.suptitle(f"Simulação CA-{model_type} | {info}", fontsize=13, fontweight='bold')

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


def plot_comparison(results_dict: dict, metric: str = 'I',
                    xlabel: str = 'Passo', ylabel: str = 'Infectados',
                    title: str = 'Comparação de Cenários',
                    save_path: str = None):
    """Compara curvas de múltiplos experimentos."""
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = plt.cm.tab10.colors

    for idx, (label, df) in enumerate(results_dict.items()):
        ax.plot(df['step'], df[metric], label=label,
                color=colors[idx % 10], linewidth=2)

    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=13)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


def plot_comparison_mc(results_mc: dict, metric_label: str = 'Infectados',
                       xlabel: str = 'Passo', title: str = 'Comparação de Cenários',
                       save_path: str = None):
    """
    Plota curvas de comparação com média e intervalo de confiança (IC 95%) de Monte Carlo.

    Args:
        results_mc: dict {label: (steps_arr, mean_arr, lo_arr, hi_arr)}
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = plt.cm.tab10.colors

    for idx, (label, (steps_arr, mean, lo, hi)) in enumerate(results_mc.items()):
        color = colors[idx % 10]
        ax.plot(steps_arr, mean, label=label, color=color, linewidth=2)
        ax.fill_between(steps_arr, lo, hi, color=color, alpha=0.18)

    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(metric_label, fontsize=11)
    ax.set_title(title, fontsize=13)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


def save_snapshots(states: list, steps: list, model_type: str,
                   output_dir: str, prefix: str = 'snapshot'):
    """Salva snapshots da grade em momentos específicos."""
    os.makedirs(output_dir, exist_ok=True)
    for state, step in zip(states, steps):
        fig = plot_grid(state, step=step, title=f'Modelo {model_type} — Passo {step}')
        path = os.path.join(output_dir, f'{prefix}_step{step:04d}.png')
        fig.savefig(path, dpi=120, bbox_inches='tight')
        plt.close(fig)
