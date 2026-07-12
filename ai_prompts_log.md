# Log de Prompts de IA — Projeto CA-Epidemiologia

Modelo utilizado: Claude Sonnet 4.6 (Anthropic)
Disciplinas: Fundamentos de Autômatos Celulares · Epidemiologia — UFRPE 2026.1

---

## 1. Brainstorm e Design do Modelo

**Prompt:**
> "Quero modelar a propagação de doenças infecciosas usando autômatos celulares.
> O modelo deve suportar SIR, SIS e SEIR em uma grade 2D toróide com vizinhanças
> de Von Neumann e Moore. Como estruturar o código em módulos Python?"

---

## 2. Correção de Código — Vetorização dos Laços

**Prompt:**
> "O código usa laços aninhados para contar vizinhos infectados, o que está lento.
> Corrija usando scipy.signal.convolve2d com boundary='wrap' para manter
> as condições de contorno toroidais."

---

## 3. Correção de Código — Análise Monte Carlo

**Prompt:**
> "Os experimentos rodam apenas uma simulação por cenário. Adicione n=30 réplicas
> Monte Carlo com sementes independentes (seed = 42 + 17r) e calcule
> média ± IC 95% para cada métrica."

---

## 4. Validação de Resultados por Log

**Prompt:**
> "Você consegue fazer logs para confirmar se o que está escrito no artigo
> está correto? Rode as simulações e compare cada valor numérico citado
> no texto com os resultados reais."

---

## 5. Correção do Artigo — Erros Técnicos

**Prompt:**
> "Analise o trabalho e traga feedbacks. Corrija: a fórmula do R0 espacial
> (usar [1-(1-β)^k]/γ no lugar de βk/γ), os valores errados da Tabela 5,
> a figura 4 ausente e o estimador de R0 retornando nan."

---

## 6. Sugestões de Melhoria Científica

**Prompt:**
> "Você poderia sugerir melhorias de qualidade científica e de escrita no artigo?
> Gostaria de tornar a linguagem mais cautelosa onde necessário, contextualizar
> melhor resultados já conhecidos na literatura, ampliar a discussão com
> interpretações e implicações práticas, e fortalecer a seção de limitações."
