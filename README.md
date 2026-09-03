# Blood Culture Contamination Tracker

> **Clinical Domain:** Diagnostic Microbiology, Antimicrobial Stewardship & Healthcare Epidemiology  
> **Standards & Guidelines:** `CLSI M47-A2`, `CLSI M100`, `CAP LAP Criteria`, `CDC NHSN BSI Guidelines`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![CLSI M47](https://img.shields.io/badge/Benchmark-Target%20%3C%203.0%25-green.svg)
![Adjudication](https://img.shields.io/badge/Algorithm-Richter%20Probability%20%2B%20DTTP-purple.svg)

</div>

---

## 📖 Clinical Overview

Blood culture contamination remains a widespread clinical challenge, leading to unnecessary broad-spectrum antimicrobial exposure (e.g., empiric vancomycin), redundant laboratory testing, repeat venipunctures, and prolonged hospital length of stay (LOS). 

The **Blood Culture Contamination Tracker** is an analytical system providing automated diagnostic decision support and infection prevention surveillance:
- **Organism Pathogenicity Adjudication:** Automated discrimination between true pathogens (e.g., *Staphylococcus aureus*, *Escherichia coli*, *Pseudomonas aeruginosa*, *Streptococcus pneumoniae*) and common skin commensals (e.g., Coagulase-negative Staphylococci, *Cutibacterium acnes*, *Corynebacterium* spp., *Bacillus* spp.).
- **Differential Time-to-Positivity (DTTP):** Evaluation of paired central venous line and peripheral blood draws for Catheter-Related Bloodstream Infection (CRBSI/CLABSI) diagnostic criteria.
- **Statistical Surveillance & Funnel Plots:** Institutional and unit-level contamination rate tracking against the **CLSI M47-A2** / **CAP** benchmark threshold of **< 3.0%**, with 95% Wilson Score confidence intervals and 2-sigma / 3-sigma funnel limits.
- **Health Economic Impact Modeling:** Quantification of avoidable direct hospital expenditure and antimicrobial overuse.

---

## 📐 Clinical Formulations & Adjudication Logic

### 1. CLSI M47-A / CAP Quality Benchmark Target
The benchmark contamination rate $R_{\text{contam}}$ is monitored against the international institutional threshold target:
$$R_{\text{contam}} = \frac{N_{\text{contaminated sets}}}{N_{\text{total sets drawn}}} \times 100\% \le 3.0\%$$

For surveillance reporting, 95% confidence intervals are estimated using the **Wilson Score Interval**:
$$CI_{95\%} = \frac{\hat{p} + \frac{z^2}{2n} \pm z \sqrt{\frac{\hat{p}(1-\hat{p})}{n} + \frac{z^2}{4n^2}}}{1 + \frac{z^2}{n}}$$
where $\hat{p} = \frac{k}{n}$, $n$ is total sets, and $z = 1.96$.

---

### 2. Richter et al. Probability Algorithm & Bottle Positivity Ratio
Probability of true bacteremia $P(\text{True Bacteremia})$ incorporates:
- **Organism Pathogenicity Weighting:** Obligate pathogens (*S. aureus*, Enterobacterales, *P. aeruginosa*, *S. pneumoniae*, *Candida* spp.) yield $>95\%$ initial prior probability. Common commensals (CoNS, diphtheroids, *Micrococcus*, *Bacillus* spp.) have low baseline probability ($\sim 5-15\%$).
- **Bottle Concordance Ratio:**
  $$\text{Ratio}_{\text{pos}} = \frac{N_{\text{pos}}}{N_{\text{total}}}$$
  - Growth in $\ge 2/2$ or $\ge 3/4$ sets/bottles of CoNS elevates the probability of true line-associated or prosthetic bacteremia ($P \approx 70-85\%$).
  - Isolated positivity in only $1/4$ or $1/2$ bottles with delayed growth strongly favors skin contamination ($P < 10\%$).

---

### 3. Differential Time-to-Positivity (DTTP) for CLABSI
When simultaneous blood cultures are obtained from an indwelling central venous catheter (CVC) and a peripheral venipuncture site:
$$\text{DTTP} = \text{TTP}_{\text{peripheral}} - \text{TTP}_{\text{central}}$$
- **CLABSI Positive Criterion:** If the central line culture becomes positive $\ge 2.0\text{ hours}$ earlier than the peripheral blood culture ($\text{DTTP} \ge 2.0\text{ h}$), catheter microbial colonization/CRBSI is verified with $>90\%$ specificity.
- **Short TTP Indicator:** Rapid growth ($\text{TTP} < 14\text{ h}$) reflects high bacterial inoculum in true sepsis, whereas prolonged growth ($\text{TTP} > 36-40\text{ h}$) with skin flora correlates with low-inoculum surface contamination.

---

### 4. Collection Site Risk Stratification
- **Peripheral Venipuncture:** Benchmark standard; sterile single-needle draw protocol.
- **Catheter / Central Line Draw:** Carries higher baseline risk of colonization and false-positive skin/hub commensal contamination unless paired peripheral cultures are drawn concurrently to establish concordance.

---

## 💻 CLI Usage & Batch Processing

The command-line interface provides tools for batch CSV surveillance, single-set clinical adjudication, and institutional surveillance statistics.

### Batch Processing Surveillance Records
Process laboratory blood culture datasets formatted according to CLSI surveillance standards:

```bash
python cli.py batch -i sample.csv -o results.csv --benchmark-target 3.0
```

#### Output Summary Example:
```text
================================================================================
  CLSI M47-A2 BLOOD CULTURE CONTAMINATION SURVEILLANCE REPORT
================================================================================
  Total Sets Evaluated:       15
  Contaminated Sets:          5
  Institutional Rate:         33.33% (Target: <= 3.0%)
  Wilson 95% Confidence Int:  [15.18%, 58.29%]
  Benchmark Status:           NON-COMPLIANT (Exceeds CLSI Benchmark)
  Estimated Excess Costs:     $35,500.00
  Potential Savings (50% red):$17,750.00
--------------------------------------------------------------------------------
  Detailed records saved to:   results.csv
================================================================================
```

### Single Set Adjudication
Adjudicate individual patient culture draws directly:

```bash
python cli.py adjudicate \
  --set-id "BC-2026-09" \
  --organism "coagulase_negative_staphylococcus" \
  --bottles-drawn 4 \
  --bottles-positive 1 \
  --ttp 46.5 \
  --site "peripheral"
```

### Surveillance Rate Calculation
Inspect an existing processed result CSV:

```bash
python cli.py surveillance -i results.csv --target 3.0
```

---

## 🐍 Python Quickstart

```python
from blood_culture_tracker import (
    BloodCultureSet,
    CultureBottle,
    CultureAdjudicationEngine,
    ContaminationSurveillanceEngine,
    EconomicImpactEngine,
)

# 1. Adjudicate paired peripheral and central line culture bottles
culture_set = BloodCultureSet(
    set_id="BC-SET-8401",
    patient_id="PT-99120",
    collection_unit="ICU",
    bottles=[
        CultureBottle("B1-CENTRAL", site_type="central_line", ttp_hours=12.5, is_positive=True, organism="coagulase_negative_staphylococcus"),
        CultureBottle("B2-PERIPH", site_type="peripheral", ttp_hours=17.0, is_positive=True, organism="coagulase_negative_staphylococcus"),
    ],
)

verdict = CultureAdjudicationEngine.adjudicate_set(culture_set)
print(f"Adjudication: {verdict['verdict']}")
print(f"Differential TTP: {verdict['dttp_hours']} hours (CLABSI suspected if >= 2.0h)")
print(f"Recommendation: {verdict['clinical_recommendation']}")

# 2. Institutional Rate & Wilson 95% CI
lower_ci, upper_ci = ContaminationSurveillanceEngine.calculate_wilson_ci(contaminated_count=5, total_sets=250)
print(f"Contamination Rate: {5/250*100:.2f}% (95% CI: [{lower_ci}%, {upper_ci}%])")

# 3. Avoidable Hospital Excess Cost Modeling
costs = EconomicImpactEngine.calculate_cost(contaminated_count=5)
print(f"Estimated Excess Cost: ${costs['estimated_annual_excess_cost_usd']:,.2f}")
```

---

## 🧪 Automated Testing

Run the full pytest suite:

```bash
python -m pytest -p no:zarr -v
```

Run the batch command verification:

```bash
python cli.py batch -i sample.csv -o out_smoke.csv
```

---

## 📄 License

MIT License. Designed for clinical microbiology laboratories, hospital infection control committees, and antimicrobial stewardship teams.

