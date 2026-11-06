#!/usr/bin/env python3
"""Recompute and verify the core BDSF Kazakhstan 2026 results.

The script deliberately reads the frozen scored dataset rather than mutable raw
exports. It writes machine-readable results and exits non-zero if a registered
result is outside its numerical tolerance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ITEMS = [f"Q{i:02d}_score" for i in range(1, 15)]
DIMENSIONS = {
    "D1": [f"Q{i:02d}_score" for i in range(1, 5)],
    "D2": [f"Q{i:02d}_score" for i in range(5, 8)],
    "D3": [f"Q{i:02d}_score" for i in range(8, 12)],
    "D4": [f"Q{i:02d}_score" for i in range(12, 15)],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def cronbach_alpha(frame: pd.DataFrame) -> float:
    values = frame.astype(float)
    k = values.shape[1]
    item_variances = values.var(axis=0, ddof=1).sum()
    total_variance = values.sum(axis=1).var(ddof=1)
    return float(k / (k - 1) * (1 - item_variances / total_variance))


def kmo_and_bartlett(frame: pd.DataFrame) -> tuple[float, float, int, float]:
    x = frame.astype(float).to_numpy()
    correlation = np.corrcoef(x, rowvar=False)
    inverse = np.linalg.inv(correlation)
    scale = np.sqrt(np.outer(np.diag(inverse), np.diag(inverse)))
    partial = -inverse / scale
    np.fill_diagonal(partial, 0.0)
    corr_off = correlation.copy()
    np.fill_diagonal(corr_off, 0.0)
    r2 = float(np.sum(corr_off**2) / 2)
    p2 = float(np.sum(partial**2) / 2)
    kmo = r2 / (r2 + p2)

    n, p = x.shape
    determinant = float(np.linalg.det(correlation))
    chi2 = -(n - 1 - (2 * p + 5) / 6) * math.log(determinant)
    df = p * (p - 1) // 2
    p_value = float(stats.chi2.sf(chi2, df))
    return kmo, chi2, df, p_value


def classify(score: pd.Series) -> pd.Series:
    rounded = score.round(2)
    return pd.cut(
        rounded,
        bins=[-np.inf, 33.32, 66.66, np.inf],
        labels=["LOW", "MEDIUM", "HIGH"],
    )


def calculate(data: pd.DataFrame) -> dict:
    included = data.loc[data["analysis_inclusion"].eq("INCLUDE")].copy()
    dimension_scores = {
        code: included[columns].mean(axis=1) for code, columns in DIMENSIONS.items()
    }
    official = pd.concat(dimension_scores, axis=1).mean(axis=1)
    levels = classify(official).value_counts().reindex(["LOW", "MEDIUM", "HIGH"], fill_value=0)
    kmo, bartlett_chi2, bartlett_df, bartlett_p = kmo_and_bartlett(included[ITEMS])

    alpha = {code: cronbach_alpha(included[columns]) for code, columns in DIMENSIONS.items()}
    alpha["items_14"] = cronbach_alpha(included[ITEMS])
    alpha["dimensions_4"] = cronbach_alpha(pd.DataFrame(dimension_scores))

    return {
        "sample": {
            "n": int(len(included)),
            "language": {str(k): int(v) for k, v in included["source_language"].value_counts().to_dict().items()},
            "gender": {str(k): int(v) for k, v in included["gender_code"].value_counts().to_dict().items()},
            "age": {str(k): int(v) for k, v in included["age_group_code"].value_counts().to_dict().items()},
        },
        "scores": {
            "BDSF_mean": float(official.mean()),
            "BDSF_median": float(official.median()),
            "BDSF_sd": float(official.std(ddof=1)),
            "D1_mean": float(dimension_scores["D1"].mean()),
            "D2_mean": float(dimension_scores["D2"].mean()),
            "D3_mean": float(dimension_scores["D3"].mean()),
            "D4_mean": float(dimension_scores["D4"].mean()),
        },
        "levels": {str(k): int(v) for k, v in levels.to_dict().items()},
        "reliability": alpha,
        "validity": {
            "KMO": kmo,
            "Bartlett_chi2": bartlett_chi2,
            "Bartlett_df": bartlett_df,
            "Bartlett_p": bartlett_p,
        },
    }


def flatten(value: object, prefix: str = "") -> dict[str, object]:
    result: dict[str, object] = {}
    if isinstance(value, dict):
        for key, child in value.items():
            child_prefix = f"{prefix}.{key}" if prefix else str(key)
            result.update(flatten(child, child_prefix))
    else:
        result[prefix] = value
    return result


def verify(actual: dict, expected: dict) -> list[dict]:
    actual_flat = flatten(actual)
    checks = []
    for metric, specification in flatten(expected).items():
        if not metric.endswith((".value", ".tolerance")):
            continue
        if not metric.endswith(".value"):
            continue
        name = metric.removesuffix(".value")
        target = specification
        tolerance = flatten(expected).get(name + ".tolerance", 0)
        observed = actual_flat.get(name)
        if isinstance(target, (int, float)) and isinstance(observed, (int, float)):
            passed = abs(float(observed) - float(target)) <= float(tolerance)
            difference = float(observed) - float(target)
        else:
            passed = observed == target
            difference = None
        checks.append({
            "metric": name,
            "observed": observed,
            "expected": target,
            "tolerance": tolerance,
            "difference": difference,
            "status": "PASS" if passed else "FAIL",
        })
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "reproduced_results")
    parser.add_argument("--expected", type=Path, default=Path(__file__).resolve().parent / "13_03_expected_results.json")
    args = parser.parse_args()

    input_path = args.project_root / "05_scoring" / "05_01_BDSF_Kazakhstan_Scored_Dataset.csv"
    if not input_path.exists():
        raise FileNotFoundError(f"Required frozen input is missing: {input_path}")
    data = pd.read_csv(input_path)
    assert len(data) == 159 and data.respondent_id.nunique() == 159
    assert data.age_group_code.isin(["18_20", "21_23", "GE27"]).all()
    assert data.analysis_inclusion.eq("INCLUDE").all()
    harmonized = pd.read_csv(args.project_root / "03_harmonization" / "03_01_BDSF_Kazakhstan_Harmonized_Dataset.csv")
    adults = harmonized.loc[harmonized.age_group_code.ne("LT18")]
    assert set(adults.respondent_id) == set(data.respondent_id)
    for q in range(1,15):
        mapped = data[f"Q{q:02d}_level"].map({"NEVER":0,"RARELY":25,"SOMETIMES":50,"OFTEN":75,"ALWAYS":100})
        assert mapped.notna().all() and np.allclose(data[f"Q{q:02d}_score"], mapped)
    actual = calculate(data)
    expected = json.loads(args.expected.read_text(encoding="utf-8"))
    checks = verify(actual, expected)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    result_document = {
        "provenance": {
            "input": str(input_path.resolve()),
            "input_sha256": sha256(input_path),
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "results": actual,
    }
    (args.output_dir / "13_07_reproduced_results.json").write_text(
        json.dumps(result_document, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    pd.DataFrame(checks).to_csv(args.output_dir / "13_08_verification_report.csv", index=False)
    failures = [row for row in checks if row["status"] == "FAIL"]
    print(f"Input SHA-256: {result_document['provenance']['input_sha256']}")
    print(f"Verification checks: {len(checks) - len(failures)}/{len(checks)} PASS")
    if failures:
        for row in failures:
            print(f"FAIL {row['metric']}: observed={row['observed']} expected={row['expected']}")
        return 1
    print("REPRODUCIBILITY STATUS: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
