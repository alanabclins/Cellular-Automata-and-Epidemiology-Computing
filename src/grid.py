"""
grid.py - Gerenciamento da grade bidimensional do autômato celular
"""

import numpy as np
from config import S, I, NEIGHBORHOOD_VON_NEUMANN, NEIGHBORHOOD_MOORE


class Grid:
    """
    Grade bidimensional para o autômato celular epidemiológico.
    Cada célula representa um indivíduo com um estado epidemiológico.
    """

    def __init__(self, size: int, neighborhood: str = NEIGHBORHOOD_MOORE, seed: int = None):
        self.size = size
        self.neighborhood = neighborhood
        if seed is not None:
            np.random.seed(seed)
        # Inicializa todos como Suscetíveis
        self.state = np.zeros((size, size), dtype=np.int8)
        self.recovery_timer = np.zeros((size, size), dtype=np.float32)
        self.exposure_timer = np.zeros((size, size), dtype=np.float32)

    def initialize(self, n_infected: int):
        """Posiciona infectados iniciais aleatoriamente na grade."""
        self.state = np.zeros((self.size, self.size), dtype=np.int8)
        self.recovery_timer = np.zeros((self.size, self.size), dtype=np.float32)
        self.exposure_timer = np.zeros((self.size, self.size), dtype=np.float32)

        positions = np.random.choice(self.size * self.size, n_infected, replace=False)
        for pos in positions:
            r, c = divmod(pos, self.size)
            self.state[r, c] = I

    def get_neighbors(self, row: int, col: int):
        """
        Retorna as coordenadas dos vizinhos de uma célula.
        Suporta condições de contorno periódicas (toroidal).
        """
        n = self.size
        if self.neighborhood == NEIGHBORHOOD_VON_NEUMANN:
            # 4 vizinhos: cima, baixo, esquerda, direita
            offsets = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        else:
            # Moore: 8 vizinhos
            offsets = [(-1, -1), (-1, 0), (-1, 1),
                       (0, -1),           (0, 1),
                       (1, -1),  (1, 0),  (1, 1)]
        return [((row + dr) % n, (col + dc) % n) for dr, dc in offsets]

    def count_infected_neighbors(self, row: int, col: int) -> int:
        """Conta quantos vizinhos estão infectados."""
        count = 0
        for r, c in self.get_neighbors(row, col):
            if self.state[r, c] == I:
                count += 1
        return count

    def get_neighbor_count(self) -> int:
        """Retorna o número de vizinhos conforme o tipo de vizinhança."""
        return 4 if self.neighborhood == NEIGHBORHOOD_VON_NEUMANN else 8

    def copy_state(self) -> np.ndarray:
        return self.state.copy()

    def get_counts(self) -> dict:
        """Retorna contagem de indivíduos em cada estado."""
        total = self.size * self.size
        unique, counts = np.unique(self.state, return_counts=True)
        result = {k: 0 for k in range(4)}
        for state, count in zip(unique, counts):
            result[int(state)] = int(count)
        return result
