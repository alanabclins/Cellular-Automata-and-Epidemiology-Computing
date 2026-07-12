# Documentação Técnica — Simulador CA-Epidemiológico

## Visão Geral

Sistema de simulação de propagação de doenças infecciosas usando Autômatos Celulares bidimensionais, implementando os modelos SIR, SIS e SEIR.

---

## Arquitetura

```
Cellular-Automata-and-Epidemiology-Computing/
├── src/
│   ├── config.py          # Constantes, estados, parâmetros padrão
│   ├── grid.py            # Grade 2D e cálculo de vizinhança
│   ├── models.py          # Regras de transição epidemiológicas
│   ├── simulation.py      # Motor de simulação
│   ├── statistics.py      # Coleta e análise de métricas
│   ├── visualization.py   # Geração de figuras
│   ├── experiments.py     # Experimentos automatizados
│   └── main.py            # Ponto de entrada CLI
├── figures/               # Figuras geradas
├── experiments/           # Resultados CSV
├── docs/                  # Documentação técnica
└── article/               # Artigo científico
```

---

## Módulos

### config.py
Define constantes globais: estados `S=0, E=1, I=2, R=3`, cores, parâmetros padrão e tipos de vizinhança.

### grid.py — Classe `Grid`

| Método | Descrição |
|--------|-----------|
| `__init__(size, neighborhood, seed)` | Inicializa grade N×N |
| `initialize(n_infected)` | Posiciona infectados iniciais aleatoriamente |
| `get_counts()` | Retorna contagem por estado |

Vizinhança toroidal: `(row + dr) % N` para condições periódicas de borda.

### models.py — Classe `EpidemiologicalModel`

Regra de infecção: `P(S→I | k vizinhos infectados) = 1 − (1 − β)^k`

| Método | Descrição |
|--------|-----------|
| `apply_sir(state, timer, grid)` | Transições SIR vetorizadas |
| `apply_sis(state, timer, grid)` | Transições SIS (R→S após recuperação) |
| `apply_seir(state, rec_t, exp_t, grid)` | Transições SEIR com incubação |
| `step(grid)` | Executa um passo de atualização síncrona |

### simulation.py — Classe `Simulation`

| Método | Descrição |
|--------|-----------|
| `run(steps, callback)` | Executa N passos com callback opcional |
| `reset(config)` | Reinicia simulação |
| `get_summary()` | Resumo estatístico |

### statistics.py — Classe `SimulationStatistics`

| Método | Descrição |
|--------|-----------|
| `peak_infection()` | Passo e contagem do pico |
| `attack_rate()` | Proporção total infectada |
| `basic_reproduction_number_estimate()` | Estimativa empírica de R₀ |
| `to_dataframe()` | Exporta histórico como DataFrame |

---

## Uso pela Linha de Comando

```bash
python main.py --mode all                          # demos + experimentos
python main.py --mode demo --model SEIR            # apenas demonstração
python main.py --mode exp                          # apenas experimentos
python main.py --mode sim --beta 0.4 --size 100   # simulação customizada
```

## Uso Programático

```python
from simulation import Simulation

cfg = {'grid_size': 100, 'beta': 0.3, 'gamma': 0.05,
       'model': 'SIR', 'neighborhood': 'moore', 'seed': 42}

sim = Simulation(cfg)
stats = sim.run()
print(stats.summary())
```

---

## Parâmetros

| Parâmetro | Intervalo | Descrição |
|-----------|-----------|-----------|
| `grid_size` | 20–500 | Dimensão N×N |
| `beta` | 0–1 | Taxa de transmissão |
| `gamma` | 0–1 | Taxa de recuperação |
| `sigma` | 0–1 | Taxa E→I (SEIR) |
| `neighborhood` | von_neumann / moore | Tipo de vizinhança |
| `model` | SIR / SIS / SEIR | Modelo epidemiológico |
| `seed` | qualquer | Semente aleatória |
