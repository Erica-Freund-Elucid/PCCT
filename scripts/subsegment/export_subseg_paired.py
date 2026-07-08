"""Export sub-segment paired data to a paired_data.csv-format file (raw totals +
Len per endpoint), reusing run_gate_analyses.load_paired_data. Lets the stat
report compute sub-segment results the same way it does canonical.

Usage: python export_subseg_paired.py <pcct_dir> <eid_dir> <out.csv>
"""
import csv, os, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
sys.path.insert(0, ROOT)
import run_gate_analyses as g  # noqa: E402


def main():
    pcct_dir, eid_dir, out = sys.argv[1], sys.argv[2], sys.argv[3]
    paired = g.load_paired_data(pcct_dir, eid_dir)
    vars_ = list(g.GATE3_PRIMARY) + list(g.GATE3_DESCRIPTIVE) + list(g.GATE3_SECONDARY)
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        hdr = ["patient_id"]
        for v in vars_:
            hdr += [f"PCCT_{v}", f"EID_{v}", f"diff_{v}", f"mean_{v}"]
        w.writerow(hdr)
        for p in paired:
            row = [p["patient_id"]]
            for v in vars_:
                pv, ev = p["pcct"].get(v), p["eid"].get(v)
                if pv is not None and ev is not None:
                    row += [pv, ev, pv - ev, (pv + ev) / 2]
                else:
                    row += ["", "", "", ""]
            w.writerow(row)
    print(f"Wrote {out}: {len(paired)} paired patients")


if __name__ == "__main__":
    main()
