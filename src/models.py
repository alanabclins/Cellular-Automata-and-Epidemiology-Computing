"""
models.py - Regras de transição para modelos epidemiológicos via Autômatos Celulares

Referência: Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals
            Capítulos 1 e 2 - Modelos compartimentais SIR e SIS
"""

import numpy as np
from scipy.signal import convolve2d
from config import S, E, I, R, MODEL_SIR, MODEL_SIS, MODEL_SEIR, NEIGHBORHOOD_VON_NEUMANN

# Kernels de vizinhança para contagem vetorizada de infectados
_KERNEL_MOORE = np.array([[1, 1, 1],
                           [1, 0, 1],
                           [1, 1, 1]], dtype=np.float32)

_KERNEL_VON_NEUMANN = np.array([[0, 1, 0],
                                 [1, 0, 1],
                                 [0, 1, 0]], dtype=np.float32)


class EpidemiologicalModel:
    """
    Classe base para modelos epidemiológicos implementados como AC.

    A lógica de transição de estados segue as equações diferenciais dos modelos
    compartimentais clássicos, adaptadas para atualização discreta e espacial:

    SIR: dS/dt = -βSI/N, dI/dt = βSI/N - γI, dR/dt = γI
    SIS: dS/dt = -βSI/N + γI, dI/dt = βSI/N - γI
    SEIR: dS/dt = -βSI/N, dE/dt = βSI/N - σE, dI/dt = σE - γI, dR/dt = γI
    """

    def __init__(self, model_type: str, beta: float, gamma: float,
                 sigma: float = 0.1, neighbor_count: int = 8):
        self.model_type = model_type
        self.beta = beta        # Taxa de transmissão por contato
        self.gamma = gamma      # Taxa de recuperação
        self.sigma = sigma      # Taxa de incubação (SEIR)
        self.neighbor_count = neighbor_count

    def infection_probability(self, n_infected_neighbors: int) -> float:
        """
        Probabilidade de uma célula suscetível se infectar dado o número de
        vizinhos infectados. Derivada da taxa β do modelo compartimental.

        P(infecção) = 1 - (1 - β)^n_infected
        """
        if n_infected_neighbors == 0:
            return 0.0
        return 1.0 - (1.0 - self.beta) ** n_infected_neighbors

    def _count_infected_neighbors(self, state: np.ndarray, neighborhood: str) -> np.ndarray:
        """
        Conta vizinhos infectados de todas as células simultaneamente via convolução 2D.
        Condições de contorno periódicas (toroidal) via boundary='wrap'.
        """
        kernel = (_KERNEL_VON_NEUMANN if neighborhood == NEIGHBORHOOD_VON_NEUMANN
                  else _KERNEL_MOORE)
        infected = (state == I).astype(np.float32)
        return convolve2d(infected, kernel, mode='same', boundary='wrap').astype(np.int32)

    def apply_sir(self, state: np.ndarray, recovery_timer: np.ndarray,
                  grid) -> tuple:
        """Aplica regras de transição do modelo SIR (vetorizado)."""
        new_state = state.copy()
        new_timer = recovery_timer.copy()
        shape = state.shape

        n_inf = self._count_infected_neighbors(state, grid.neighborhood)

        # S → I: probabilidade local P = 1 - (1-β)^k
        s_mask = (state == S)
        p_inf = np.where(n_inf > 0, 1.0 - (1.0 - self.beta) ** n_inf, 0.0)
        new_infections = s_mask & (np.random.random(shape) < p_inf)
        new_state[new_infections] = I
        new_timer[new_infections] = 0.0

        # I → R
        i_mask = (state == I)
        new_timer[i_mask] += 1.0
        new_recoveries = i_mask & (np.random.random(shape) < self.gamma)
        new_state[new_recoveries] = R
        new_timer[new_recoveries] = 0.0  # zera timer ao sair do estado

        return new_state, new_timer

    def apply_sis(self, state: np.ndarray, recovery_timer: np.ndarray,
                  grid) -> tuple:
        """
        Aplica regras de transição do modelo SIS (vetorizado).
        Indivíduos recuperados voltam a ser suscetíveis.
        """
        new_state = state.copy()
        new_timer = recovery_timer.copy()
        shape = state.shape

        n_inf = self._count_infected_neighbors(state, grid.neighborhood)

        # S → I
        s_mask = (state == S)
        p_inf = np.where(n_inf > 0, 1.0 - (1.0 - self.beta) ** n_inf, 0.0)
        new_infections = s_mask & (np.random.random(shape) < p_inf)
        new_state[new_infections] = I
        new_timer[new_infections] = 0.0

        # I → S (sem imunidade)
        i_mask = (state == I)
        new_timer[i_mask] += 1.0
        new_recoveries = i_mask & (np.random.random(shape) < self.gamma)
        new_state[new_recoveries] = S
        new_timer[new_recoveries] = 0.0

        return new_state, new_timer

    def apply_seir(self, state: np.ndarray, recovery_timer: np.ndarray,
                   exposure_timer: np.ndarray, grid) -> tuple:
        """
        Aplica regras do modelo SEIR com período de incubação (E = Exposto) (vetorizado).
        """
        new_state = state.copy()
        new_recovery = recovery_timer.copy()
        new_exposure = exposure_timer.copy()
        shape = state.shape

        n_inf = self._count_infected_neighbors(state, grid.neighborhood)

        # S → E
        s_mask = (state == S)
        p_inf = np.where(n_inf > 0, 1.0 - (1.0 - self.beta) ** n_inf, 0.0)
        new_exposures = s_mask & (np.random.random(shape) < p_inf)
        new_state[new_exposures] = E
        new_exposure[new_exposures] = 0.0

        # E → I
        e_mask = (state == E)
        new_exposure[e_mask] += 1.0
        new_infections = e_mask & (np.random.random(shape) < self.sigma)
        new_state[new_infections] = I
        new_recovery[new_infections] = 0.0
        new_exposure[new_infections] = 0.0

        # I → R
        i_mask = (state == I)
        new_recovery[i_mask] += 1.0
        new_recoveries = i_mask & (np.random.random(shape) < self.gamma)
        new_state[new_recoveries] = R
        new_recovery[new_recoveries] = 0.0

        return new_state, new_recovery, new_exposure

    def step(self, grid) -> dict:
        """Executa um passo de atualização síncrona do autômato celular."""
        if self.model_type == MODEL_SIR:
            new_state, new_timer = self.apply_sir(
                grid.state, grid.recovery_timer, grid)
            grid.state = new_state
            grid.recovery_timer = new_timer

        elif self.model_type == MODEL_SIS:
            new_state, new_timer = self.apply_sis(
                grid.state, grid.recovery_timer, grid)
            grid.state = new_state
            grid.recovery_timer = new_timer

        elif self.model_type == MODEL_SEIR:
            new_state, new_recovery, new_exposure = self.apply_seir(
                grid.state, grid.recovery_timer, grid.exposure_timer, grid)
            grid.state = new_state
            grid.recovery_timer = new_recovery
            grid.exposure_timer = new_exposure

        return grid.get_counts()
