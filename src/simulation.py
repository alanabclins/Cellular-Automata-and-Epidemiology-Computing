"""
simulation.py - Motor principal da simulação do autômato celular epidemiológico
"""

import numpy as np
from grid import Grid
from models import EpidemiologicalModel
from statistics import SimulationStatistics
from config import DEFAULT_CONFIG


class Simulation:
    """
    Orquestra a simulação completa do modelo CA-epidemiológico.
    
    Integra a grade (Grid), o modelo (EpidemiologicalModel) e as estatísticas,
    implementando o ciclo de atualização síncrona característico dos autômatos celulares.
    """

    def __init__(self, config: dict = None):
        self.config = {**DEFAULT_CONFIG, **(config or {})}
        self.grid = None
        self.model = None
        self.stats = None
        self.running = False
        self.paused = False
        self._step = 0
        self._initialize()

    def _initialize(self):
        cfg = self.config
        self.grid = Grid(
            size=cfg['grid_size'],
            neighborhood=cfg['neighborhood'],
            seed=cfg['seed']
        )
        self.grid.initialize(n_infected=cfg['initial_infected'])

        self.model = EpidemiologicalModel(
            model_type=cfg['model'],
            beta=cfg['beta'],
            gamma=cfg['gamma'],
            sigma=cfg['sigma'],
            neighbor_count=self.grid.get_neighbor_count()
        )

        self.stats = SimulationStatistics(
            grid_size=cfg['grid_size'],
            model_type=cfg['model'],
            gamma=cfg.get('gamma'),
        )
        # Registra estado inicial
        self.stats.record(self.grid.get_counts())
        self._step = 0
        self.running = False
        self.paused = False

    def reset(self, config: dict = None):
        if config:
            self.config = {**self.config, **config}
        self._initialize()

    def step(self) -> dict:
        """Executa um único passo de simulação."""
        counts = self.model.step(self.grid)
        self.stats.record(counts)
        self._step += 1
        return counts

    def run(self, steps: int = None, callback=None) -> SimulationStatistics:
        """
        Executa a simulação por 'steps' iterações.
        callback(step, counts, grid_state) é chamado a cada passo se fornecido.
        """
        max_steps = steps or self.config['steps']
        self.running = True

        for _ in range(max_steps):
            if not self.running:
                break
            counts = self.step()
            if callback:
                callback(self._step, counts, self.grid.state.copy())
            # Para cedo se a epidemia se extinguiu
            if counts.get(2, 0) == 0 and self._step > 1:
                break

        self.running = False
        return self.stats

    def stop(self):
        self.running = False

    def pause(self):
        self.paused = True

    def resume(self):
        self.paused = False

    @property
    def current_step(self) -> int:
        return self._step

    @property
    def current_state(self) -> np.ndarray:
        return self.grid.state.copy()

    def get_summary(self) -> dict:
        return self.stats.summary()


