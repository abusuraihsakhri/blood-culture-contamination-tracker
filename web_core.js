(function (root, factory) {
  const api = factory();
  if (typeof module === "object" && module.exports) module.exports = api;
  if (root) root.BCTracker = api;
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  "use strict";

  const TRUE_PATHOGENS = new Set([
    "staphylococcus_aureus", "escherichia_coli", "klebsiella_pneumoniae",
    "pseudomonas_aeruginosa", "streptococcus_pneumoniae", "enterococcus_faecalis",
    "enterococcus_faecium", "candida_albicans", "candida_glabrata",
    "streptococcus_pyogenes", "streptococcus_agalactiae", "neisseria_meningitidis",
    "bacteroides_fragilis", "listeria_monocytogenes", "salmonella_enterica",
    "acinetobacter_baumannii", "serratia_marcescens", "proteus_mirabilis"
  ]);

  const COMMON_COMMENSALS = new Set([
    "coagulase_negative_staphylococcus", "staphylococcus_epidermidis",
    "staphylococcus_hominis", "staphylococcus_capitis", "staphylococcus_warneri",
    "corynebacterium_species", "cutibacterium_acnes", "propionibacterium_acnes",
    "micrococcus_luteus", "bacillus_species", "viridans_group_streptococci",
    "aerococcus_viridans"
  ]);

  const ALIASES = new Map([
    ["coagulase_negative_staphylococci", "coagulase_negative_staphylococcus"],
    ["cons", "coagulase_negative_staphylococcus"],
    ["corynebacterium_spp", "corynebacterium_species"],
    ["bacillus_spp", "bacillus_species"],
    ["viridans_streptococci", "viridans_group_streptococci"],
    ["propionibacterium_acnes", "cutibacterium_acnes"]
  ]);

  function normalizeOrganism(value) {
    const text = String(value ?? "")
      .trim()
      .toLowerCase()
      .replace(/[^a-z0-9]+/g, "_")
      .replace(/^_+|_+$/g, "");
    return ALIASES.get(text) || text;
  }

  function finiteNumber(value, name) {
    const number = Number(value);
    if (!Number.isFinite(number)) throw new Error(`${name} must be a finite number.`);
    return number;
  }

  function optionalNumber(value, name) {
    if (value === "" || value === null || value === undefined) return null;
    const number = finiteNumber(value, name);
    if (number < 0) throw new Error(`${name} must be non-negative.`);
    return number;
  }

  function adjudicateCase(input) {
    const organism = normalizeOrganism(input.organism);
    if (!organism) throw new Error("Enter an organism.");

    const drawn = finiteNumber(input.bottlesDrawn, "Bottles drawn");
    const positive = finiteNumber(input.bottlesPositive, "Bottles positive");
    if (!Number.isInteger(drawn) || drawn <= 0) throw new Error("Bottles drawn must be a positive integer.");
    if (!Number.isInteger(positive) || positive < 0 || positive > drawn) {
      throw new Error("Bottles positive must be an integer between 0 and bottles drawn.");
    }

    const ttp = optionalNumber(input.ttpHours, "TTP");
    const peripheral = optionalNumber(input.peripheralTtpHours, "Peripheral TTP");
    const central = optionalNumber(input.centralTtpHours, "Central-line TTP");
    if (peripheral !== null && central !== null && positive < 2) {
      throw new Error("Paired central/peripheral TTP values require at least two positive bottles/cultures.");
    }

    if (positive === 0) {
      return {
        organism,
        verdict: "Negative culture / no positive bottles represented",
        isContamination: false,
        score: 0,
        concordance: 0,
        dttpHours: null,
        reasons: ["No positive culture bottle is represented."],
        recommendation: "Interpret with the complete clinical and microbiology context."
      };
    }

    const isPathogen = TRUE_PATHOGENS.has(organism);
    const isCommensal = COMMON_COMMENSALS.has(organism);
    let score = 0;
    const reasons = [];

    if (isPathogen) {
      score += 3;
      reasons.push("Organism is in the configured pathogen list (+3).");
    } else if (isCommensal) {
      score -= 2;
      reasons.push("Organism is in the configured common-commensal list (-2).");
    } else {
      reasons.push("Organism is not classified by the built-in organism lists.");
    }

    const shortest = [ttp, peripheral, central].filter((value) => value !== null);
    const minTtp = shortest.length ? Math.min(...shortest) : null;
    if (minTtp !== null && minTtp < 14) {
      score += 2;
      reasons.push(`Rapid time to positivity (${minTtp.toFixed(1)} h < 14 h) increases the true-BSI signal (+2).`);
    } else if (minTtp !== null && minTtp > 36 && isCommensal) {
      score -= 2;
      reasons.push(`Delayed time to positivity (${minTtp.toFixed(1)} h > 36 h) with a common commensal decreases the true-BSI signal (-2).`);
    }

    const concordance = positive / drawn;
    if (concordance >= 0.75 && drawn >= 2) {
      score += 2;
      reasons.push(`High bottle concordance (${positive}/${drawn} positive) (+2).`);
    } else if (concordance <= 0.25 && drawn >= 4) {
      score -= 2;
      reasons.push(`Low bottle concordance (${positive}/${drawn} positive) (-2).`);
    }

    let dttpHours = null;
    let catheterSource = false;
    if (central !== null && peripheral !== null) {
      dttpHours = Math.round((peripheral - central) * 100) / 100;
      if (dttpHours >= 2) {
        catheterSource = true;
        score += 3;
        reasons.push(`Central-line culture was positive ${dttpHours.toFixed(2)} h earlier; DTTP ≥ 2 h supports a catheter source (+3).`);
      } else if (dttpHours <= -2) {
        reasons.push(`Peripheral culture was positive ${Math.abs(dttpHours).toFixed(2)} h earlier; DTTP does not support a catheter source.`);
      } else {
        reasons.push(`Paired DTTP was ${dttpHours.toFixed(2)} h, below the 2 h catheter-source threshold.`);
      }
    }

    let verdict;
    let isContamination = false;
    let recommendation;
    if (catheterSource) {
      verdict = "Catheter-related bloodstream infection (CRBSI) supported by DTTP";
      recommendation = "Correlate clinically and apply local CRBSI guidance. This is not an NHSN CLABSI determination.";
    } else if (score >= 2) {
      verdict = "Likely true bloodstream infection signal";
      recommendation = "The rule set favors true bloodstream infection; correlate with clinical findings and laboratory context.";
    } else if (score <= -2) {
      verdict = "Probable blood-culture contamination signal";
      isContamination = true;
      recommendation = "The rule set favors contamination; review repeat cultures, collection details, and clinical context before changing therapy.";
    } else {
      verdict = "Indeterminate; clinical correlation and repeat culture may be needed";
      recommendation = "Use clinical assessment, repeat cultures when appropriate, and local policy.";
    }

    return { organism, verdict, isContamination, score, concordance, dttpHours, reasons, recommendation };
  }

  function wilsonCI(contaminated, total, z = 1.96) {
    const k = finiteNumber(contaminated, "Contaminated count");
    const n = finiteNumber(total, "Total count");
    if (!Number.isInteger(k) || !Number.isInteger(n) || n < 0 || k < 0 || k > n) {
      throw new Error("Counts must be integers satisfying 0 ≤ contaminated ≤ total.");
    }
    if (n === 0) return [0, 0];
    if (!(z > 0)) throw new Error("z must be positive.");
    const p = k / n;
    const z2 = z * z;
    const denominator = 1 + z2 / n;
    const center = (p + z2 / (2 * n)) / denominator;
    const margin = z * Math.sqrt((p * (1 - p) + z2 / (4 * n)) / n) / denominator;
    return [Math.max(0, center - margin) * 100, Math.min(1, center + margin) * 100];
  }

  function surveillance(input) {
    const total = finiteNumber(input.total, "Total cultures");
    const contaminated = finiteNumber(input.contaminated, "Contaminated cultures");
    const target = finiteNumber(input.targetPct, "Target");
    if (!Number.isInteger(total) || !Number.isInteger(contaminated) || total < 0 || contaminated < 0 || contaminated > total) {
      throw new Error("Counts must be integers satisfying 0 ≤ contaminated ≤ total.");
    }
    if (target < 0 || target > 100) throw new Error("Target must be between 0 and 100%.");
    const rate = total ? (contaminated / total) * 100 : 0;
    const ci = wilsonCI(contaminated, total);
    return { total, contaminated, ratePct: rate, ci95Pct: ci, targetPct: target, meetsTarget: rate <= target };
  }

  return Object.freeze({
    TRUE_PATHOGENS,
    COMMON_COMMENSALS,
    normalizeOrganism,
    adjudicateCase,
    wilsonCI,
    surveillance
  });
});
