"""
statistics.py - Coleta e análise de dados epidemiológicos da simulação
"""

import numpy as np
import pandas as pd
from config import S, E, I, R


class SimulationStatistics:
    """Coleta e armazena estatísticas ao longo da simulação."""

    def __init__(self, grid_size: int, model_type: str, gamma: float = None):
        self.grid_size = grid_size
        self.total_population = grid_size * grid_size
        self.model_type = model_type
        self.gamma = gamma
        self.history = []
        self.step_count = 0

    def record(self, counts: dict):
        """Registra contagens de um passo de simulação."""
        record = {
            'step': self.step_count,
            'S': counts.get(S, 0),
            'E': counts.get(E, 0),
            'I': counts.get(I, 0),
            'R': counts.get(R, 0),
        }
        record['prevalence'] = record['I'] / self.total_population
        self.history.append(record)
        self.step_count += 1

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.history)

    def peak_infection(self) -> dict:
        df = self.to_dataframe()
        idx = df['I'].idxmax()
        return {'step': int(df.loc[idx, 'step']), 'count': int(df.loc[idx, 'I'])}

    def epidemic_duration(self) -> int:
        """Retorna o passo em que a infecção se extingue."""
        df = self.to_dataframe()
        extinct = df[df['I'] == 0]
        if len(extinct) > 0:
            return int(extinct.iloc[0]['step'])
        return self.step_count

    def total_infected(self) -> int:
        """Total acumulado de infectados (= R final no SIR)."""
        df = self.to_dataframe()
        if self.model_type == 'SIR':
            return int(df['R'].iloc[-1])
        else:
            return int(df['I'].max() + df.get('R', pd.Series([0])).iloc[-1])

    def attack_rate(self) -> float:
        """Proporção da população que foi infectada."""
        return self.total_infected() / self.total_population

    def basic_reproduction_number_estimate(self) -> float:
        """
        Estimativa empírica do R₀ pela taxa de crescimento exponencial inicial.

        Na fase inicial (I pequeno, S ≈ N), I(t) ≈ I(0)·exp(r·t) onde
        r = β·k - γ é a taxa líquida de crescimento. A relação R₀ = 1 + r/γ
        (aproximação de campo médio para a fase exponencial) permite estimar R₀
        mesmo quando S_final = 0, caso em que a relação de tamanho final diverge.

        Usa regressão log-linear sobre os primeiros passos de crescimento (até o pico).
        """
        df = self.to_dataframe()
        peak_idx = df['I'].idxmax()
        # Usa no mínimo 5 e no máximo 20 pontos da fase de crescimento
        growth = df.iloc[1:max(peak_idx, 6)]
        growth = growth[growth['I'] > 0]
        if len(growth) < 3:
            return float('nan')
        try:
            t = growth['step'].values.astype(float)
            log_i = np.log(growth['I'].values.astype(float))
            # Regressão linear: log I(t) = log I(0) + r·t
            coeffs = np.polyfit(t, log_i, 1)
            r = coeffs[0]            # taxa de crescimento log por passo
            rho = np.exp(r)          # fator de crescimento por passo
            if self.gamma and self.gamma > 0:
                # Relação exata para AC discreto: ρ = 1 + γ·(R₀_AC - 1)
                # → R₀_AC = 1 + (ρ - 1) / γ
                r0 = 1.0 + (rho - 1.0) / self.gamma
            else:
                r0 = rho             # aproximação se γ não disponível
            return round(max(0.0, float(r0)), 3)
        except Exception:
            return float('nan')

    def summary(self) -> dict:
        return {
            'model': self.model_type,
            'total_steps': self.step_count,
            'peak_infection': self.peak_infection(),
            'epidemic_duration': self.epidemic_duration(),
            'total_infected': self.total_infected(),
            'attack_rate': round(self.attack_rate(), 4),
            'estimated_R0': self.basic_reproduction_number_estimate(),
        }
