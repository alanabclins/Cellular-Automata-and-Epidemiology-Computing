# Documento de Requisitos — Simulador CA-Epidemiológico

**Versão:** 1.0 · **UFRPE 2026.1**

---

## Requisitos Funcionais

- **RF01** — Configuração: grade, β, γ, vizinhança, modelo, infectados iniciais
- **RF02** — Controle: executar, pausar, reiniciar, múltiplas simulações
- **RF03** — Modelos: SIR, SIS, SEIR com regras probabilísticas locais
- **RF04** — AC: atualização síncrona, Von Neumann, Moore, toroidal
- **RF05** — Visualização: mapa 2D, curvas epidemiológicas, snapshots, PNG
- **RF06** — Estatísticas: pico, duração, taxa de ataque, R₀, CSV
- **RF07** — Experimentos: variação de β, vizinhança, modelo, grade, I₀

---

## Requisitos Não Funcionais

- **RNF01** — Desempenho: grade 100×100, 200 passos em menos de 60s
- **RNF02** — Reprodutibilidade: semente aleatória fixa
- **RNF03** — Modularidade: módulos independentes
- **RNF04** — Documentação: docstrings, README, docs/
- **RNF05** — Tecnologia: Python 3.8+, NumPy, Matplotlib, Pandas, SciPy

---

## Rastreabilidade

| Requisito | Módulo | Status |
|-----------|--------|--------|
| RF01 | config.py, main.py | ✅ |
| RF02 | simulation.py | ✅ |
| RF03 | models.py | ✅ |
| RF04 | grid.py | ✅ |
| RF05 | visualization.py | ✅ |
| RF06 | statistics.py | ✅ |
| RF07 | experiments.py | ✅ |
| RNF01–05 | Todos | ✅ |
