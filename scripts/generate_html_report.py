"""Generate a self-contained, editable HTML report of PCCT qualification gate findings.

Reads gate_results/, tracker/, and the workflow + case-summary spreadsheets.
Writes gate_results/qualification_report.html (self-contained, all images embedded).

The HTML has text marked contenteditable="true" so it can be edited live in
the browser. A "Download edited HTML" button serializes the current DOM
back to a file. Numeric data blocks (gate_summary.txt content) are kept
non-editable.

Run from PCCT/ root or from PCCT/scripts/:
    python scripts/generate_html_report.py
"""
import base64
import csv
import io
import os
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "gate_results" / "qualification_report.html"


def b64_image(p: Path) -> str:
    if not p.exists():
        return ""
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode("ascii")


def fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


def read(p: Path) -> str:
    return p.read_text(encoding="utf-8") if p.exists() else ""


def extract_block(text: str, header_re: str, end_re: str = r"^={5,}") -> str:
    m = re.search(header_re, text, re.MULTILINE)
    if not m:
        return ""
    rest = text[m.end():]
    end = re.search(end_re, rest, re.MULTILINE)
    return rest[: end.start()].strip() if end else rest.strip()


def html_escape(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


# ── Plot builders ─────────────────────────────────────────────────────────────

def load_snr_csv(p: Path):
    """Return dict pid -> {mean_hu, std, snr} (skip rows with err)."""
    out = {}
    for r in csv.DictReader(open(p, encoding="utf-8")):
        if r.get("err"):
            continue
        try:
            pid = r["case"].split("_Bv_")[0]
            out[pid] = {
                "mean_hu": float(r["mean"]),
                "std": float(r["std"]),
                "snr": float(r["snr"]),
            }
        except (ValueError, KeyError):
            pass
    return out


def plot_contrast_timing(pcct, eid, threshold=250):
    """Paired-bar plot: per-patient aortic mean HU for PCCT and EID, threshold line."""
    paired_pids = sorted(set(pcct) & set(eid))
    pcct_vals = [pcct[p]["mean_hu"] for p in paired_pids]
    eid_vals = [eid[p]["mean_hu"] for p in paired_pids]

    fig, ax = plt.subplots(figsize=(11, 4.5))
    x = np.arange(len(paired_pids))
    w = 0.38
    ax.bar(x - w/2, pcct_vals, w, label=f"PCCT (N={len(pcct_vals)})", color="#2563eb", edgecolor="white")
    ax.bar(x + w/2, eid_vals, w, label=f"EID (N={len(eid_vals)})", color="#f59e0b", edgecolor="white")
    ax.axhline(threshold, color="#dc2626", linestyle="--", linewidth=1.2,
               label=f"Threshold {threshold} HU")
    ax.set_xticks(x)
    ax.set_xticklabels(paired_pids, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Mean aortic HU (10 mm ROI at centroid)")
    ax.set_title("Gate 1.2 Contrast Timing — paired PCCT vs EID aortic enhancement")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    ax.set_ylim(bottom=0)
    return fig


def plot_snr_paired(pcct, eid):
    """Side-by-side: paired bar chart of SNR, and scatter."""
    paired_pids = sorted(set(pcct) & set(eid))
    p = [pcct[i]["snr"] for i in paired_pids]
    e = [eid[i]["snr"] for i in paired_pids]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.5), gridspec_kw={"width_ratios": [2.4, 1]})

    # Left: paired bars
    x = np.arange(len(paired_pids)); w = 0.38
    ax1.bar(x - w/2, p, w, label="PCCT", color="#2563eb", edgecolor="white")
    ax1.bar(x + w/2, e, w, label="EID", color="#f59e0b", edgecolor="white")
    ax1.set_xticks(x); ax1.set_xticklabels(paired_pids, rotation=60, ha="right", fontsize=8)
    ax1.set_ylabel("SNR (mean / std HU)")
    ax1.set_title("Gate 1.3 Image Noise — paired SNR (aortic ROI)")
    ax1.legend(fontsize=9); ax1.grid(axis="y", alpha=0.25)

    # Right: scatter with y=x
    mx = max(max(p), max(e)) * 1.1
    ax2.plot([0, mx], [0, mx], "k--", alpha=0.4, linewidth=1, label="y=x (parity)")
    ax2.scatter(e, p, c="#2563eb", s=44, zorder=5, edgecolor="white")
    ax2.set_xlabel("EID SNR")
    ax2.set_ylabel("PCCT SNR")
    ax2.set_title(f"PCCT vs EID SNR\nratio {np.mean([a/b for a,b in zip(p,e)]):.2f}× (PCCT/EID)")
    ax2.set_xlim(0, mx); ax2.set_ylim(0, mx); ax2.set_aspect("equal")
    ax2.legend(fontsize=8); ax2.grid(alpha=0.25)

    fig.tight_layout()
    return fig


def load_effort():
    """Return (paired_pids, pcct_eff, ccta_eff) for patients with numeric effort on both."""
    import openpyxl
    wb = openpyxl.load_workbook(ROOT / "PCCT_CCTA_Case_Summaries.xlsx", data_only=True)
    ws = wb["Case Summaries"]
    hdr = [ws.cell(2, c).value for c in range(1, ws.max_column + 1)]
    eff_col = hdr.index("Analyst Effort\n(1=Most, 5=Least)") + 1

    by_pid = defaultdict(dict)
    for r in range(3, ws.max_row + 1):
        pid = ws.cell(r, 1).value
        scan = ws.cell(r, 2).value
        eff = ws.cell(r, eff_col).value
        if not pid or scan not in ("PCCT", "CCTA"):
            continue
        try:
            e = float(eff)
        except (TypeError, ValueError):
            continue
        by_pid[pid][scan] = e

    pids = sorted([pid for pid, d in by_pid.items() if "PCCT" in d and "CCTA" in d])
    return pids, [by_pid[p]["PCCT"] for p in pids], [by_pid[p]["CCTA"] for p in pids]


def plot_effort(pids, pcct_eff, ccta_eff):
    """Two-panel: paired effort bars, and histogram of deltas."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.5), gridspec_kw={"width_ratios": [2.4, 1]})

    x = np.arange(len(pids)); w = 0.38
    ax1.bar(x - w/2, pcct_eff, w, label="PCCT", color="#2563eb", edgecolor="white")
    ax1.bar(x + w/2, ccta_eff, w, label="CCTA", color="#f59e0b", edgecolor="white")
    ax1.set_xticks(x); ax1.set_xticklabels(pids, rotation=60, ha="right", fontsize=8)
    ax1.set_ylabel("Analyst effort (1=most, 5=least)")
    ax1.set_title(f"Gate 2.4 Analyst Editing Effort — paired (N={len(pids)})")
    ax1.set_yticks([1, 2, 3, 4, 5])
    ax1.set_ylim(0, 5.5)
    ax1.legend(fontsize=9); ax1.grid(axis="y", alpha=0.25)

    deltas = [p - c for p, c in zip(pcct_eff, ccta_eff)]
    mean_d = np.mean(deltas)
    bins = np.arange(-4.25, 4.25, 0.5)
    ax2.hist(deltas, bins=bins, color="#2563eb", edgecolor="white", alpha=0.8)
    ax2.axvline(0, color="black", linewidth=1, alpha=0.5)
    ax2.axvline(mean_d, color="#dc2626", linewidth=1.5, label=f"mean {mean_d:+.2f}")
    ax2.set_xlabel("Effort delta (PCCT − CCTA)\n← PCCT harder    PCCT easier →")
    ax2.set_ylabel("Patients")
    ax2.set_title(f"Effort delta distribution\nmedian {np.median(deltas):+.1f}")
    ax2.legend(fontsize=9); ax2.grid(alpha=0.25)

    fig.tight_layout()
    return fig


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    gate_summary = read(ROOT / "gate_results" / "gate_summary.txt")
    gate1 = extract_block(gate_summary, r"^GATE 1.*?$")
    gate2 = extract_block(gate_summary, r"^GATE 2.*?$")
    gate3 = extract_block(gate_summary, r"^GATE 3.*?$")
    gate4 = extract_block(gate_summary, r"^GATE 4 .{,2}SYSTEMATIC", end_re=r"^GATE 4 SUPPLEMENTARY")

    # BA plots (canonical only)
    plot_dir = ROOT / "gate_results" / "bland_altman_plots"
    ba_vars = ["LumenVol", "WallVol", "VesselVol", "CALCVol", "LRNCVol",
               "NonCALCMATXVol", "TotalPlaqueVolume"]
    ba_canonical = {v: b64_image(plot_dir / f"BA_{v}.png") for v in ba_vars}

    # New Gate 1 / Gate 2.4 plots — generate from raw data
    snr_pcct = load_snr_csv(ROOT / "gate_results" / "snr_pcct.csv")
    snr_eid = load_snr_csv(ROOT / "gate_results" / "snr_eid.csv")
    contrast_img = fig_to_b64(plot_contrast_timing(snr_pcct, snr_eid))
    snr_img = fig_to_b64(plot_snr_paired(snr_pcct, snr_eid))

    eff_pids, eff_pcct, eff_ccta = load_effort()
    effort_img = fig_to_b64(plot_effort(eff_pids, eff_pcct, eff_ccta))

    # Audit + case reviews
    audit_md = read(ROOT / "tracker" / "workitem-audit.md")
    m = re.search(r"## Summary\s*\n\n(.*?)\n\n", audit_md, re.DOTALL)
    audit_summary_md = m.group(1) if m else ""
    case_reviews_md = read(ROOT / "tracker" / "case-reviews.md")

    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>PCCT Qualification - Gate Findings Report</title>
<style>
:root {{
    --bg: #fafbfc; --panel: #fff; --border: #d1d5db; --text: #1f2937;
    --muted: #6b7280; --pass: #16a34a; --fail: #dc2626; --warn: #f59e0b; --accent: #2563eb;
}}
body {{ background: var(--bg); color: var(--text); font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; line-height: 1.55; max-width: 1100px; margin: 0 auto; padding: 24px; padding-top: 64px; }}
h1, h2, h3, h4 {{ color: var(--text); margin-top: 1.4em; }}
h1 {{ font-size: 1.85em; border-bottom: 2px solid var(--accent); padding-bottom: 8px; }}
h2 {{ font-size: 1.4em; border-bottom: 1px solid var(--border); padding-bottom: 4px; }}
h3 {{ font-size: 1.15em; color: var(--accent); }}
.meta {{ color: var(--muted); font-size: 0.9em; }}
.panel {{ background: var(--panel); border: 1px solid var(--border); border-radius: 6px; padding: 16px 20px; margin: 12px 0; }}
.pass {{ color: var(--pass); font-weight: 600; }}
.fail {{ color: var(--fail); font-weight: 600; }}
.warn {{ color: var(--warn); font-weight: 600; }}
table {{ border-collapse: collapse; width: 100%; margin: 8px 0 16px; font-size: 0.92em; }}
th {{ background: #f3f4f6; text-align: left; padding: 8px 10px; border-bottom: 2px solid var(--border); }}
td {{ padding: 6px 10px; border-bottom: 1px solid #e5e7eb; vertical-align: top; }}
tr:hover {{ background: #f9fafb; }}
pre {{ background: #f3f4f6; padding: 12px; border-radius: 4px; overflow-x: auto; font-size: 0.85em; font-family: Consolas, Menlo, monospace; white-space: pre-wrap; }}
.callout {{ border-left: 4px solid var(--accent); background: #eff6ff; padding: 10px 14px; margin: 10px 0; }}
.callout-warn {{ border-left-color: var(--warn); background: #fffbeb; }}
img.fig {{ max-width: 100%; height: auto; border: 1px solid var(--border); border-radius: 4px; margin: 6px 0; }}
.ba-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(360px, 1fr)); gap: 12px; }}
.ba-grid figure {{ margin: 0; }}
.ba-grid figcaption {{ font-size: 0.85em; color: var(--muted); margin-top: 4px; }}
#editbar {{ position: fixed; top: 0; left: 0; right: 0; background: #1f2937; color: white; padding: 10px 20px; z-index: 1000; box-shadow: 0 2px 4px rgba(0,0,0,0.1); display: flex; align-items: center; gap: 12px; }}
#editbar button {{ background: var(--accent); color: white; border: 0; padding: 6px 14px; border-radius: 4px; cursor: pointer; font-size: 0.92em; }}
#editbar button:hover {{ background: #1d4ed8; }}
#editbar .info {{ font-size: 0.85em; opacity: 0.85; }}
[contenteditable="true"]:focus {{ outline: 2px solid var(--accent); outline-offset: 2px; background: #fefce8; }}
[contenteditable="true"]:hover {{ background: #fef9c3; }}
.uneditable {{ background: #f3f4f6; padding: 8px 12px; border-radius: 4px; font-size: 0.85em; color: var(--muted); border-left: 3px solid var(--muted); }}
</style>
</head>
<body>

<div id="editbar">
    <strong>PCCT Qualification Report</strong>
    <span class="info">Yellow-highlighted text is editable — click to edit</span>
    <span style="flex:1"></span>
    <button onclick="downloadEdited()">Download edited HTML</button>
    <button onclick="window.print()">Print / PDF</button>
</div>

<h1 contenteditable="true">PCCT Qualification — Gate Findings Report</h1>
<p class="meta" contenteditable="true">Generated {now} from <code>gate_results/gate_summary.txt</code>, <code>tracker/*.md</code>, paired SNR data, and analyst case-summary spreadsheet. Reference: B.1P Delta Validation OQ (730-CVV-040 v0.1) and Original B.1P OQ (4-B1P-033 v2.0).</p>

<section>
<h2>Executive Summary</h2>
<div class="panel" contenteditable="true">
<p><strong>N = 25 paired patients</strong> (vessel-overlap restricted; PT-124 excluded due to non-overlapping vessel territories). Preliminary results — target N ≥ 30.</p>
<table>
<tr><th>Gate</th><th>Requirement</th><th>Status</th></tr>
<tr><td><strong>Gate 1</strong> — Technical & Image Quality</td><td>DICOM compliance, contrast timing, SNR, kernel</td><td><span class="pass">PASS</span></td></tr>
<tr><td><strong>Gate 2</strong> — Workflow Integration</td><td>Ingestion, centerline, lumen/wall editing</td><td><span class="warn">REVIEW</span> (centerline coverage diff)</td></tr>
<tr><td><strong>Gate 3</strong> — Quantitative Reproducibility</td><td>95% CI overlap with delta OQ wCV</td><td><span class="pass">PASS</span> (all 7 endpoints)</td></tr>
<tr><td><strong>Gate 4</strong> — Bias & Agreement</td><td>Bland-Altman bias + proportional bias</td><td><span class="warn">MIXED</span> — process outputs PASS, plaque components FAIL on bias</td></tr>
</table>
</div>
</section>

<section>
<h2>Methodology Notes</h2>
<div class="callout" contenteditable="true">
<strong>Vessel-overlap normalization (canonical as of 2026-05-08).</strong> For each PCCT/EID pair, summaries are summed only over the (bodySite, vessel-location) intersection — so the same anatomical extent is measured on each scanner. Patients with no overlapping vessel are excluded (PT-124). Most pairs are partial-overlap; PCCT typically traces 1-3 additional distal vessels (PDA, PLB, marginals, diagonals).
</div>
<div class="callout" contenteditable="true">
<strong>Metric per endpoint family.</strong> Log-wCV (and log-scale Bland-Altman) for process outputs — Lumen, Wall, Vessel Volume — matching 4-B1P-033 Table 7 and 730-CVV-040 Table 6. Untransformed wCV (and untransformed BA) for plaque volumes (CALC, LRNC, NonCALC Matrix, Total Plaque), matching 4-B1P-033 Table 9 reporting scale.
</div>
<div class="callout callout-warn" contenteditable="true">
<strong>Wall bias threshold is project-specific.</strong> |bias| &lt; 10% of mean was added 2026-05-08 in response to case-review finding (PT-142) of likely true modality bias for Wall at high disease. Not derivable from 4-B1P-033 or 730-CVV-040 (both assess Wall reproducibility by wCV only). Threshold matched to plaque-component bias thresholds since Wall segmentation is dominated by the wall-plaque boundary.
</div>
</section>

<section>
<h2>Gate 1 — Technical & Image Quality Prerequisites</h2>

<h3>1.2 Contrast Timing — paired aortic enhancement</h3>
<img class="fig" src="{contrast_img}" alt="Contrast timing paired bars">
<p class="meta" contenteditable="true">Threshold: peak aortic HU ≥ 250 in ≥ 90% of cases. PCCT 28/28 (100%) PASS. EID 24/25 (96%) PASS. Only PT-124 (242 HU) below threshold.</p>

<h3>1.3 Image Noise — paired SNR</h3>
<img class="fig" src="{snr_img}" alt="SNR paired and scatter">
<p class="meta" contenteditable="true">SNR = mean_HU / std_HU within a 10 mm aortic ROI at centroid. Threshold: PCCT SNR non-inferior to EID (ratio ≥ 0.85). PCCT mean SNR ~2× EID; lower noise in 24/25 paired cases. PT-142 and PT-158 EID excluded — Aorta.nrrd not generated.</p>

<pre class="uneditable">{html_escape(gate1)}</pre>
<div class="panel" contenteditable="true">
<strong>Reviewer commentary:</strong> Editable space for Gate 1 notes.
</div>
</section>

<section>
<h2>Gate 2 — Workflow Integration</h2>
<pre class="uneditable">{html_escape(gate2)}</pre>

<h3>2.4 Lumen & Wall Editing — analyst effort comparison</h3>
<img class="fig" src="{effort_img}" alt="Effort paired bars and delta histogram">
<p class="meta" contenteditable="true">Per-case analyst effort (1=most effort, 5=least), recorded in PCCT_CCTA_Case_Summaries.xlsx. Paired N={len(eff_pids)} patients with numeric scores on both scans. PCCT mean {np.mean(eff_pcct):.2f}, CCTA mean {np.mean(eff_ccta):.2f}; delta +{np.mean(eff_pcct)-np.mean(eff_ccta):.2f} (median {np.median([p-c for p,c in zip(eff_pcct, eff_ccta)]):.1f}). Symmetric distribution → PCCT does not impose meaningfully higher editing effort.</p>

<div class="panel" contenteditable="true">
<strong>Reviewer commentary:</strong> Centerline coverage and editing-effort findings.
</div>
</section>

<section>
<h2>Gate 3 — Quantitative Reproducibility</h2>
<pre class="uneditable">{html_escape(gate3)}</pre>
<div class="panel" contenteditable="true">
<strong>Reviewer commentary:</strong> All 7 endpoints PASS at N=25 with vessel-overlap normalization. Wall recovered from FAIL on target-overlap to PASS on vessel-overlap, confirming the prior fail was a length-normalization artifact from PCCT's longer distal trace covering thin-walled regions.
</div>
</section>

<section>
<h2>Gate 4 — Bias & Agreement</h2>
<pre class="uneditable">{html_escape(gate4)}</pre>

<h3>Bland-Altman plots</h3>
<div class="ba-grid">
"""
    for v in ba_vars:
        if ba_canonical[v]:
            html += f'<figure><img class="fig" src="{ba_canonical[v]}" alt="BA {v}"><figcaption contenteditable="true">{v}</figcaption></figure>\n'
    html += """</div>

<div class="panel" contenteditable="true">
<strong>Reviewer commentary:</strong> Lumen, Wall, Vessel, Total Plaque bias all PASS. CALC, LRNC, NonCALC component biases FAIL — directional pattern consistent with case-review notes (PT-142 high-disease Wall modality bias).
</div>
</section>

<section>
<h2>Case Reviews</h2>
<div class="panel" contenteditable="true">
"""
    cr = case_reviews_md
    cr = re.sub(r"^# (.*)$", r"<h3>\1</h3>", cr, flags=re.MULTILINE)
    cr = re.sub(r"^### (.*)$", r"<h4>\1</h4>", cr, flags=re.MULTILINE)
    cr = re.sub(r"^---$", r"<hr>", cr, flags=re.MULTILINE)
    cr = re.sub(r"^\- \*\*(.+?):\*\* (.*)$", r"<p><strong>\1:</strong> \2</p>", cr, flags=re.MULTILINE)
    cr = re.sub(r"^\- (.*)$", r"<li>\1</li>", cr, flags=re.MULTILINE)
    cr = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", cr)
    cr = re.sub(r"`(.+?)`", r"<code>\1</code>", cr)
    html += cr
    html += """</div>
</section>

<section>
<h2>Workitem Audit Summary</h2>
<div class="panel" contenteditable="true">
<p>From <code>tracker/workitem-audit.md</code> — audits all 81 workitems in the workflow spreadsheet against the ip3-manager instance and local <code>workitem_summaries/</code> state.</p>
"""
    audit_lines = audit_summary_md.split("\n")
    if audit_lines:
        html += "<table>\n"
        for line in audit_lines:
            if not line.strip() or line.strip().startswith("|---"):
                continue
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            tag = "th" if cells[0] == "Status" else "td"
            html += "<tr>" + "".join(f"<{tag}>{c}</{tag}>" for c in cells) + "</tr>\n"
        html += "</table>\n"
    html += """<p><em>Carolyn Walker's 13 MISSING_LOCAL workitems are intentionally not in the gate-analysis cohort.</em></p>
</div>
</section>

<script>
function downloadEdited() {
    const clone = document.documentElement.cloneNode(true);
    const bar = clone.querySelector('#editbar');
    if (bar) bar.remove();
    const html = '<!DOCTYPE html>\\n' + clone.outerHTML;
    const blob = new Blob([html], {type: 'text/html'});
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    const ts = new Date().toISOString().replace(/[:.]/g, '-').slice(0,16);
    a.download = `qualification_report_edited_${ts}.html`;
    a.click();
    URL.revokeObjectURL(url);
}
</script>
</body>
</html>
"""

    OUT.write_text(html, encoding="utf-8")
    print(f"Wrote {OUT}")
    print(f"  size: {len(html.encode('utf-8'))//1024} KB")


if __name__ == "__main__":
    main()
