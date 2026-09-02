# Blood Culture Contamination Tracker

> **Domain:** Infectious Disease Surveillance & Microbiology  
> **Reference Guidelines & Standards:** `CLSI M100, EUCAST & CDC NHSN Clinical Standards`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

Blood Culture Contamination Tracker
CLSI contamination rate (target <3%) with organism-specific true vs contaminant adjudication.
Stdlib parser / mapper with batch CSV and single lookup.

Blood Culture Contamination Tracker & TTP Differential Engine
=============================================================
Clinical microbiology surveillance and diagnostic decision support adhering to
CLSI M47-A2 and CDC NHSN bloodstream infection standards:
- Time-To-Positivity (TTP) and Differential TTP (DTTP) for CLABSI diagnosis
- Organism pathogenicity matrix & bottle concordance adjudication
- Wilson score confidence intervals and Statistical Process Control (SPC) funnel limits
- Hospital economic impact & avoidable excess cost modeling

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`AdjudicationVerdict`** — dedicated module for adjudication verdict evaluation and state verification.
- **`CultureBottle`** — dedicated module for culture bottle evaluation and state verification.
- **`BloodCultureSet`** — dedicated module for blood culture set evaluation and state verification.
- **`CultureAdjudicationEngine`**: Classifies blood culture positivity as true bacteremia, contaminant, or CLABSI.
- **`ContaminationSurveillanceEngine`**: Computes contamination rates, Wilson confidence intervals, and unit-level funnel limits.
- **`EconomicImpactEngine`**: Calculates excess hospital costs and avoidable antimicrobial expenditures.

---

## 📐 Mathematical Formulation & Logic

```text
  score = 0
  score = 0  # Positive = true bacteremia signal, Negative = contamination signal
  Calculate 95% Wilson Score Confidence Interval for a proportion:
  ci_lower, ci_upper = cls.calculate_wilson_ci(contaminated_sets, total_sets)
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --task-id <value> --target <value> --primary <value> --secondary <value>
```

### Parameter Reference
- `--task-id`: Specifies input measurement or parameter value.
- `--target`: Specifies input measurement or parameter value.
- `--primary`: Specifies input measurement or parameter value.
- `--secondary`: Specifies input measurement or parameter value.
- `--critical`: Specifies input measurement or parameter value.
- `--status`: Specifies input measurement or parameter value.
- `--input`: Specifies input measurement or parameter value.
- `--output`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `query` | Parameter / observation metric | Required |
| `name` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t blood-culture-contamination-tracker .
docker run -p 8000:8000 blood-culture-contamination-tracker
```
