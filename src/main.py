"""
main.py - Ponto de entrada principal do simulador CA-Epidemiológico
Uso: python main.py [--mode sim|exp|demo] [--model SIR|SIS|SEIR]
"""

import sys
import os
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from simulation import Simulation
from visualization import plot_combined, save_snapshots
from experiments import run_all_experiments
from config import DEFAULT_CONFIG, MODEL_SIR, MODEL_SIS, MODEL_SEIR


def run_demo(model: str = MODEL_SIR, output_dir: str = '../figures'):
    """Executa uma simulação demonstração e salva figuras."""
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n=== Simulação Demonstração — Modelo {model} ===")

    cfg = {
        **DEFAULT_CONFIG,
        'model': model,
        'grid_size': 80,
        'beta': 0.3,
        'gamma': 0.05,
        'steps': 300,
        'seed': 42,
        'initial_infected': 5,
    }

    sim = Simulation(cfg)
    snap_states, snap_steps = [], []
    N = cfg['grid_size'] ** 2

    def collect(step, counts, state):
        # Snapshots em marcos específicos
        milestones = {10, 30, 60, 100, 150, 200}
        if step in milestones:
            snap_states.append(state.copy())
            snap_steps.append(step)

    print("Executando simulação...")
    stats = sim.run(callback=collect)
    df = stats.to_dataframe()
    summary = stats.summary()

    print(f"\nResultados:")
    print(f"  Modelo:              {summary['model']}")
    print(f"  Pico de Infectados:  {summary['peak_infection']['count']} "
          f"({summary['peak_infection']['count']/N*100:.1f}%) "
          f"no passo {summary['peak_infection']['step']}")
    print(f"  Duração da Epidemia: {summary['epidemic_duration']} passos")
    print(f"  Total Infectados:    {summary['total_infected']} "
          f"({summary['attack_rate']*100:.1f}%)")
    print(f"  R₀ estimado:         {summary['estimated_R0']}")

    # Figura combinada final
    final_path = os.path.join(output_dir, f'demo_{model}_final.png')
    fig = plot_combined(
        state=sim.current_state,
        df=df,
        step=sim.current_step,
        model_type=model,
        config=cfg,
        save_path=final_path
    )
    plt.close(fig)
    print(f"\nFigura final salva em: {os.path.abspath(final_path)}")

    # Snapshots
    save_snapshots(snap_states, snap_steps, model,
                   output_dir=output_dir, prefix=f'snap_{model}')
    print(f"Snapshots salvos em: {os.path.abspath(output_dir)}")

    return summary


def main():
    parser = argparse.ArgumentParser(
        description='Simulador CA-Epidemiológico — Autômatos Celulares + Epidemiologia')
    parser.add_argument('--mode', choices=['sim', 'exp', 'demo', 'all'],
                        default='all', help='Modo de execução')
    parser.add_argument('--model', choices=['SIR', 'SIS', 'SEIR'],
                        default='SIR', help='Modelo epidemiológico')
    parser.add_argument('--beta', type=float, default=None)
    parser.add_argument('--gamma', type=float, default=None)
    parser.add_argument('--size', type=int, default=None)
    parser.add_argument('--steps', type=int, default=None)
    parser.add_argument('--output', type=str, default='../figures')

    args = parser.parse_args()

    output_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), args.output)
    os.makedirs(output_dir, exist_ok=True)

    if args.mode in ('demo', 'all'):
        for model in [MODEL_SIR, MODEL_SIS, MODEL_SEIR]:
            run_demo(model=model, output_dir=output_dir)

    if args.mode in ('exp', 'all'):
        tables = run_all_experiments()
        # Salva tabelas CSV
        exp_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), '../experiments')
        os.makedirs(exp_dir, exist_ok=True)
        for name, df in tables.items():
            df.to_csv(os.path.join(exp_dir, f'results_{name}.csv'), index=False)
        print(f"\nTabelas CSV salvas em: {os.path.abspath(exp_dir)}")

    if args.mode == 'sim':
        cfg = {**DEFAULT_CONFIG, 'model': args.model}
        if args.beta: cfg['beta'] = args.beta
        if args.gamma: cfg['gamma'] = args.gamma
        if args.size: cfg['grid_size'] = args.size
        if args.steps: cfg['steps'] = args.steps
        run_demo(model=args.model, output_dir=output_dir)


if __name__ == '__main__':
    main()
