"""
config.py - Configurações globais do simulador CA-Epidemiológico
"""

# Estados epidemiológicos
S = 0   # Suscetível
E = 1   # Exposto (para SEIR)
I = 2   # Infectado
R = 3   # Recuperado

# Cores para visualização (RGB normalizado 0-1)
STATE_COLORS = {
    S: (0.2, 0.7, 0.2),    # Verde - Suscetível
    E: (1.0, 0.8, 0.0),    # Amarelo - Exposto
    I: (0.85, 0.1, 0.1),   # Vermelho - Infectado
    R: (0.1, 0.3, 0.85),   # Azul - Recuperado
}

STATE_LABELS = {
    S: 'Suscetível',
    E: 'Exposto',
    I: 'Infectado',
    R: 'Recuperado',
}

# Tipos de vizinhança
NEIGHBORHOOD_VON_NEUMANN = 'von_neumann'
NEIGHBORHOOD_MOORE = 'moore'

# Modelos epidemiológicos
MODEL_SIR = 'SIR'
MODEL_SIS = 'SIS'
MODEL_SEIR = 'SEIR'

# Parâmetros padrão
DEFAULT_CONFIG = {
    'grid_size': 100,
    'beta': 0.3,            # Taxa de transmissão
    'gamma': 0.05,          # Taxa de recuperação
    'sigma': 0.1,           # Taxa de incubação (SEIR)
    'initial_infected': 5,
    'neighborhood': NEIGHBORHOOD_MOORE,
    'model': MODEL_SIR,
    'steps': 200,
    'seed': 42,
}
