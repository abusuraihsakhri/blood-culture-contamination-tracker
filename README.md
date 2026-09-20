# Blood Culture Contamination Tracker

### [Open the Live Application →](https://abusuraihsakhri.github.io/blood-culture-contamination-tracker/)

A small Python and browser-based utility for transparent blood-culture contamination heuristics, differential time-to-positivity (DTTP) interpretation, Wilson confidence intervals, and contamination-rate surveillance.

The calculations are rule-based decision-support aids. They are not a validated diagnostic model and do not replace laboratory, infectious-disease, infection-prevention, or NHSN surveillance review.

## Features

- Classifies represented organisms using a small, explicit pathogen/common-commensal list.
- Combines organism category, time to positivity, and bottle concordance into an inspectable rule score.
- Calculates paired DTTP as `peripheral TTP - central-line TTP`; a matching-organism central culture becoming positive at least 2 hours earlier supports a catheter source/CRBSI assessment.
- Calculates overall contamination rate and Wilson 95% confidence intervals.
- Provides a configurable contamination-rate quality target; 3% is the default for continuity with the historically cited benchmark, not a universal pass/fail standard.
- Includes an optional assumption-based cost calculator. The bundled dollar values are illustrative defaults and should be replaced with local estimates.
- Runs as a zero-runtime-dependency Python CLI or as a static browser application. The browser application uses JavaScript directly; Pyodide is not required.

## Browser application

Open `index.html` locally, or use the deployed GitHub Pages site once enabled for the repository. The interface is light by default, includes a dark-mode toggle, and runs entirely in the browser.

The browser UI does not request patient names, medical-record numbers, dates of birth, or other identifiers. It does not send case inputs to a server. Only the selected theme is stored in browser `localStorage`.

## CLI

Python 3.10 or newer is recommended. Runtime code uses only the standard library.

Batch-process the sample CSV:

```bash
python cli.py batch -i sample.csv -o results.csv
```

Review one represented culture set:

```bash
python cli.py adjudicate \
  --organism coagulase_negative_staphylococcus \
  --bottles-drawn 2 \
  --bottles-positive 2 \
  --central-ttp 14 \
  --peripheral-ttp 18 \
  --site central_line
```

Calculate surveillance statistics from a processed CSV:

```bash
python cli.py surveillance -i results.csv --target 3.0
```

Run the deterministic stress test:

```bash
python simulator.py 500 --seed 7
```

## Input notes

Organism names are normalized to lowercase underscore-separated identifiers. The built-in lists are intentionally limited and should not be treated as a complete microbiology taxonomy. Unknown organisms remain unclassified rather than being matched by substring.

For paired DTTP, central-line and peripheral cultures should be collected in a clinically appropriate paired context and yield the same organism. DTTP is a clinical CRBSI/catheter-source concept; it is not, by itself, a CDC/NHSN CLABSI surveillance criterion.

The contamination-rate target is configurable because local definitions, collection practices, and quality goals differ. A historically cited overall contamination threshold of 3% remains common, while lower rates such as 1% are achievable in high-performing programs.

## Testing

Install the test dependency and run:

```bash
python -m pip install pytest
python -m pytest -p no:zarr -v
python -m compileall -q .
node --check web_core.js
node --check app.js
```

GitHub Actions tests Python 3.10, 3.11, and 3.12, exercises the CLI and simulator, checks browser JavaScript, and deploys the static site from `master` after tests pass.

## Project structure

- `blood_culture_tracker.py` — canonical adjudication, surveillance, Wilson CI, and cost-model calculations.
- `cli.py` — validated single-case, batch CSV, and surveillance commands.
- `ttp_differential.py` — compatibility helpers delegated to the canonical core.
- `blood_culture.py` — small organism-list lookup compatibility CLI.
- `simulator.py` — deterministic stress test of the real adjudication path.
- `web_core.js`, `app.js`, `index.html`, `styles.css` — static GitHub Pages application.
- `tests/` — regression tests.

## References and terminology

- CDC, *Background Information: Terminology and Estimates of Risk* — distinguishes clinical CRBSI from NHSN CLABSI surveillance terminology: https://www.cdc.gov/infection-control/hcp/intravascular-catheter-related-infection/terminology-estimates-of-risk.html
- CDC NHSN, *Bloodstream Infection Event (Central Line-Associated Bloodstream Infection and Non-central Line Associated Bloodstream Infection)*: https://www.cdc.gov/nhsn/psc/bsi/index.html
- IDSA/SCCM guideline, *Evaluating New Fever in Adult Patients in the ICU* — discusses paired blood cultures and differential time to positivity for catheter-associated bacteremia: https://www.idsociety.org/practice-guideline/new-fever-in-critically-ill-patients/
- CDC, *Blood Culture Contamination: An Overview for Infection Control and Antibiotic Stewardship Programs Working with the Clinical Laboratory*: https://www.cdc.gov/antibiotic-use/media/pdfs/fs-bloodculture-508.pdf

## Privacy and intended use

The repository contains no external API integration and requires no secrets. Browser calculations are local. Do not place protected health information in public repositories, issue reports, or shared example files. Validate the rules and thresholds against current local policy and clinical guidance before operational use.

## License

MIT License. See [LICENSE](LICENSE).
