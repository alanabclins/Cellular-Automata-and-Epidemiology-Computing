# Cellular Automata and Epidemiology Computing

Simulador epidemiológico baseado em autômatos celulares (AC) implementando os modelos SIR, SIS e SEIR em uma grade toroide 2D. Desenvolvido como trabalho final das disciplinas de **Fundamentos de Autômatos Celulares** e **Epidemiologia** — UFRPE, 2026.1.

**Autora:** Alana Barbalho Camara Lins  
**Orientador:** Prof. Jones Oliveira de Albuquerque  

---

## Estrutura do repositório

```
.
├── src/                  # Código-fonte do simulador
│   ├── main.py           # Ponto de entrada
│   ├── simulation.py     # Motor principal da simulação
│   ├── models.py         # Regras de transição SIR / SIS / SEIR
│   ├── grid.py           # Grade toroide e vizinhanças
│   ├── experiments.py    # Experimentos Monte Carlo
│   ├── statistics.py     # Coleta e análise de métricas
│   ├── visualization.py  # Geração de gráficos
│   └── config.py         # Parâmetros padrão
│
├── experiments/          # Resultados em CSV (n = 30 réplicas por cenário)
│   ├── results_beta.csv
│   ├── results_neighborhood.csv
│   ├── results_models.csv
│   ├── results_population.csv
│   └── results_initial.csv
│
├── figures/              # Figuras geradas pelos experimentos
│
└── article/
    ├── Modelagem_Espacial_CA_Epidemiologia_UFRPE2026.pdf  # Artigo em PDF
    └── latex/            # Fonte LaTeX (Nature Scientific Reports)
        ├── main.tex
        ├── sample.bib
        └── wlscirep.cls  (+ arquivos de estilo)
```

---

## Como executar

```bash
# 1. Clone o repositório
git clone https://github.com/alanabclins/Cellular-Automata-and-Epidemiology-Computing
cd Cellular-Automata-and-Epidemiology-Computing/src

# 2. Instale as dependências
pip install numpy scipy matplotlib pandas

# 3. Execute todos os experimentos
python3 main.py --mode all

# 4. Execute apenas uma simulação interativa
python3 main.py --mode demo
```

As figuras serão salvas em `figures/` e os CSVs em `experiments/`.

---

## Modelo

Cada célula da grade representa um indivíduo em um dos estados `{S, E, I, R}`. A atualização é síncrona e a probabilidade de infecção depende do número de vizinhos infectados:

$$P(S \to I) = 1 - (1 - \beta)^{k_I}$$

O número de reprodução espacial sob as hipóteses do modelo é:

$$R_0^{AC} = \frac{1-(1-\beta)^k}{\gamma}$$

que difere do valor de campo médio $R_0^{ODE} = \beta/\gamma$ e satura em $1/\gamma$ quando $\beta \to 1$.

---

## Parâmetros de referência (influenza sazonal)

| Parâmetro | Valor | Descrição |
|-----------|-------|-----------|
| β | 0,18 | Taxa de transmissão por vizinho por passo |
| γ | 0,14 | Probabilidade de recuperação por passo (~7 dias) |
| σ | 0,50 | Probabilidade de progressão E→I (~2 dias) |
| N | 80×80 | Tamanho da grade |
| k | 8 (Moore) | Número de vizinhos |

---

## Artigo

O manuscrito completo está disponível em dois formatos:

- **PDF:** [`article/Modelagem_Espacial_CA_Epidemiologia_UFRPE2026.pdf`](article/Modelagem_Espacial_CA_Epidemiologia_UFRPE2026.pdf)
- **LaTeX:** pasta `article/latex/`, compilável diretamente no [Overleaf](https://www.overleaf.com) via upload da pasta

---

## Uso de IA

O desenvolvimento do simulador contou com auxílio do modelo Claude Sonnet 4.6 (Anthropic). O histórico de prompts utilizado está disponível neste repositório.
