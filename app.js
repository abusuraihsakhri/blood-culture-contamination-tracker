"use strict";

document.addEventListener("DOMContentLoaded", () => {
  const core = window.BCTracker;
  const root = document.documentElement;
  const themeButton = document.getElementById("theme-toggle");

  function setTheme(theme) {
    root.dataset.theme = theme;
    localStorage.setItem("bc-tracker-theme", theme);
    themeButton.setAttribute("aria-label", theme === "dark" ? "Switch to light theme" : "Switch to dark theme");
  }

  setTheme(localStorage.getItem("bc-tracker-theme") === "dark" ? "dark" : "light");
  themeButton.addEventListener("click", () => setTheme(root.dataset.theme === "dark" ? "light" : "dark"));

  document.querySelectorAll(".tab").forEach((button) => {
    button.addEventListener("click", () => {
      document.querySelectorAll(".tab").forEach((tab) => {
        const active = tab === button;
        tab.classList.toggle("is-active", active);
        tab.setAttribute("aria-selected", String(active));
      });
      const tabName = button.dataset.tab;
      document.getElementById("case-panel").classList.toggle("is-hidden", tabName !== "case");
      document.getElementById("surveillance-panel").classList.toggle("is-hidden", tabName !== "surveillance");
    });
  });

  function value(id) { return document.getElementById(id).value; }
  function setText(id, text) { document.getElementById(id).textContent = text; }

  function setChip(id, text, state) {
    const chip = document.getElementById(id);
    chip.textContent = text;
    chip.className = `status-chip ${state}`;
  }

  document.getElementById("case-example").addEventListener("click", () => {
    document.getElementById("organism").value = "coagulase_negative_staphylococcus";
    document.getElementById("bottles-drawn").value = "2";
    document.getElementById("bottles-positive").value = "2";
    document.getElementById("ttp").value = "14";
    document.getElementById("draw-site").value = "central_line";
    document.getElementById("central-ttp").value = "14";
    document.getElementById("peripheral-ttp").value = "18";
    document.getElementById("case-form").requestSubmit();
  });

  document.getElementById("case-form").addEventListener("submit", (event) => {
    event.preventDefault();
    setText("case-error", "");
    try {
      const result = core.adjudicateCase({
        organism: value("organism"),
        bottlesDrawn: value("bottles-drawn"),
        bottlesPositive: value("bottles-positive"),
        ttpHours: value("ttp"),
        centralTtpHours: value("central-ttp"),
        peripheralTtpHours: value("peripheral-ttp")
      });
      setText("case-verdict", result.verdict);
      setText("case-score", String(result.score));
      setText("case-concordance", `${Math.round(result.concordance * 100)}%`);
      setText("case-dttp", result.dttpHours === null ? "—" : `${result.dttpHours.toFixed(2)} h`);
      setText("case-recommendation", result.recommendation);
      setChip(
        "contamination-chip",
        result.isContamination ? "Contamination signal" : "Not flagged as contamination",
        result.isContamination ? "warn" : "good"
      );

      const list = document.getElementById("case-reasons");
      list.replaceChildren();
      result.reasons.forEach((reason) => {
        const item = document.createElement("li");
        item.textContent = reason;
        list.appendChild(item);
      });
    } catch (error) {
      setText("case-error", error.message || "Unable to analyze the supplied values.");
    }
  });

  document.getElementById("surveillance-form").addEventListener("submit", (event) => {
    event.preventDefault();
    setText("surveillance-error", "");
    try {
      const result = core.surveillance({
        total: value("total-cultures"),
        contaminated: value("contaminated-cultures"),
        targetPct: value("target-pct")
      });
      const [lower, upper] = result.ci95Pct;
      setText("surveillance-rate", `${result.ratePct.toFixed(2)}%`);
      setText("surveillance-ci", `${lower.toFixed(2)}–${upper.toFixed(2)}%`);
      setText("surveillance-target", `≤ ${result.targetPct.toFixed(1)}%`);
      setText("surveillance-status", result.meetsTarget ? "Within configured target" : "Above configured target");
      setChip("target-chip", result.meetsTarget ? "Within target" : "Above target", result.meetsTarget ? "good" : "warn");
      setText(
        "surveillance-copy",
        `${result.contaminated} of ${result.total} sets were classified as contaminated. The Wilson 95% confidence interval is ${lower.toFixed(2)}% to ${upper.toFixed(2)}%.`
      );
    } catch (error) {
      setText("surveillance-error", error.message || "Unable to analyze the supplied values.");
    }
  });
});
