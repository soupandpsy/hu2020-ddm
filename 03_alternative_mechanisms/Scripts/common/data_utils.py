"""Audited Matching choice and four-condition coding; no data exclusion."""
import numpy as np
import pandas as pd
from common.baseline_v import CONDITION_DESIGN, DESIGN_COLUMNS


def encode_conditions(data):
    """Return a fresh frame with a unique index and the shared design columns.

    Existing event IDs and source rows survive index reset. Values are validated,
    never silently repaired or converted from accuracy coding.
    """
    required = {"val", "id", "stim", "response"}
    if not required.issubset(data.columns):
        raise ValueError("Missing coding columns: " + str(sorted(required-set(data.columns))))
    if data[list(required)].isna().any().any():
        raise ValueError("Missing condition, stimulus or choice values")
    if not set(data.stim.unique()).issubset({0, 1}) or not set(data.response.unique()).issubset({0, 1}):
        raise ValueError("stim and response must be numeric 0/1, not accuracy labels")
    pairs = list(zip(data.val, data.id))
    if any(pair not in CONDITION_DESIGN for pair in pairs):
        raise ValueError("Expected Good/Bad crossed with Self/Other")
    result = data.copy().reset_index(drop=True)
    matrix = np.asarray([CONDITION_DESIGN[pair] for pair in pairs], dtype=float).reshape(-1, len(DESIGN_COLUMNS))
    for j, name in enumerate(DESIGN_COLUMNS):
        if name in result and not np.array_equal(result[name].to_numpy(), matrix[:, j]):
            raise ValueError("Existing design column disagrees: " + name)
        result[name] = matrix[:, j]
    sign = 2 * result.stim.to_numpy() - 1
    if "sign" in result and not np.array_equal(result.sign.to_numpy(), sign):
        raise ValueError("Existing sign disagrees with stimulus")
    result["sign"] = sign
    result["match"] = result.stim.map({0:"mismatch",1:"match"})
    return result


def validate_choice_mapping(data):
    """Check saved choice, accuracy and subject-specific key mapping when present."""
    if "accuracy" in data:
        if not data.accuracy.isin([0, 1]).all() or not np.array_equal(
                (data.response == data.stim).astype(int), data.accuracy):
            raise ValueError("Accuracy and choice/stimulus coding disagree")
    if "match" in data:
        if not data["match"].isin(["match", "mismatch"]).all() or not np.array_equal(
                (data["match"] == "match").astype(int), data.stim):
            raise ValueError("Match label and stimulus disagree")
    if {"response_key", "upper_key"}.issubset(data.columns):
        if not data.response_key.isin(["m", "n"]).all() or not data.upper_key.isin(["m", "n"]).all():
            raise ValueError("Non-task key in valid-response data")
        if not np.array_equal((data.response_key == data.upper_key).astype(int), data.response):
            raise ValueError("Choice and actual key mapping disagree")
        if not data.groupby("subj_idx").upper_key.nunique().eq(1).all():
            raise ValueError("Inconsistent upper key within subject")
    return True


"""Valid-response subset I/O. This is NOT a deadline/censoring observation model."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
REPOSITORY = PROJECT.parent
SCOPE = "valid_response_subset_unadjusted_for_unknown_deadline"
RT_MIN_SECONDS = 0.2
RT_POLICY = "rt_ge_0p2_20261107"


def validate_data(data, require_design=True):
    required = ["subj_idx", "rt", "response", "stim", "val", "id"]
    missing = set(required)-set(data.columns)
    if missing:
        raise ValueError("Missing columns: " + str(sorted(missing)))
    if data.empty or data[required].isna().any().any():
        raise ValueError("Data must be nonempty with no missing required values")
    if not pd.api.types.is_numeric_dtype(data.rt) or not np.isfinite(data.rt).all() or not data.rt.gt(0).all():
        raise ValueError("RT must be finite positive seconds")
    if not pd.api.types.is_numeric_dtype(data.subj_idx) or not np.isfinite(data.subj_idx).all() or not (data.subj_idx % 1 == 0).all():
        raise ValueError("Subject identifiers must be finite integers")
    encode_conditions(data)
    validate_choice_mapping(data)
    for col in ("event_id", "source_row"):
        if col in data and (data[col].isna().any() or data[col].duplicated().any()):
            raise ValueError(col + " must be nonmissing and globally unique")
    if require_design:
        expected = {(v, i, s) for v in ("Good", "Bad") for i in ("Self", "Other") for s in (0, 1)}
        for subject, frame in data.groupby("subj_idx"):
            if set(zip(frame.val, frame.id, frame.stim)) != expected:
                raise ValueError("Incomplete four-condition x two-direction design for " + str(subject))
    return True


def prepare_data(data, require_design=True):
    """Prepare an already selected valid-response frame; never drop records.

    Raw experimental tables are deliberately rejected: selection is inherited
    through the audited saved input, not reimplemented as a new exclusion rule.
    """
    validate_data(data, require_design=require_design)
    result = encode_conditions(data)
    if "event_id" not in result:
        if "source_row" in result:
            result["event_id"] = "source:" + result.source_row.astype(str)
        else:
            result["event_id"] = ["row:" + str(i) for i in range(len(result))]
    result.attrs["analysis_scope"] = SCOPE
    return result


def dataset_hash(data):
    """Hash exact ordered trial identity, values and columns, excluding attrs/index."""
    frame = data.copy().reset_index(drop=True)
    payload = frame.to_csv(index=False, float_format="%.17g", lineterminator="\n")
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    def convert(value):
        if isinstance(value, Path): return str(value)
        if isinstance(value, np.generic): return value.item()
        if isinstance(value, np.ndarray): return value.tolist()
        raise TypeError("Not JSON serializable: " + type(value).__name__)
    text = json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False, default=convert)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text + "\n", encoding="utf-8")
    temporary.replace(path)


def load_data(path=None, *, rt_min=RT_MIN_SECONDS):
    """Load CSV. Default checks existing audited originals and reconstructs lineage.

    Empirical input excludes RT < 0.2 seconds, for both response accuracies.
    Pass rt_min=None explicitly when reading historical or simulated data.
    No timeout removal, unit guessing or fresh subject exclusion is performed.
    """
    default = REPOSITORY / "02_reanalysis_2026/Data/matching_input.csv"
    path = default if path is None else Path(path)
    frame = pd.read_csv(path, float_precision="round_trip")
    if path.resolve() == default.resolve():
        audit_path = REPOSITORY / "02_reanalysis_2026/Output/data_audit.json"
        audit = json.loads(audit_path.read_text(encoding="utf-8"))["matching"]
        original = REPOSITORY / "01_reproduction_2020/Data/matching_input.csv"
        raw_path = REPOSITORY / "01_reproduction_2020/Data/raw/confirm/MS_rep_matchingTask_raw.csv"
        for f, key in ((original, "original_sha256"), (raw_path, "raw_sha256")):
            if hashlib.sha256(f.read_bytes()).hexdigest() != audit[key]:
                raise ValueError("Audited source hash changed: " + str(f))
        if len(frame) != audit["fit_rows"] or frame.subj_idx.nunique() != audit["subjects"]:
            raise ValueError("Audited row/subject counts changed")
        raw = pd.read_csv(raw_path)
        if "source_row" not in frame or not frame.source_row.between(2, len(raw)+1).all():
            raise ValueError("Missing or invalid audited source-row lineage")
        selected = raw.iloc[frame.source_row.to_numpy(dtype=int)-2].reset_index(drop=True)
        if not selected.CorrResp.isin(["m", "n"]).all():
            raise ValueError("Unexpected correct response key in audited source")
        upper = selected.CorrResp.where(selected.Match == "match", selected.CorrResp.map({"m": "n", "n": "m"}))
        expected = {
            "subj_idx": selected.Subject, "val": selected.Morality,
            "id": selected.Identity, "rt": selected.RT,
            "stim": (selected.Match == "match").astype(int),
            "response_key": selected.Resp, "upper_key": upper,
            "accuracy": selected.ACC, "response": (selected.Resp == upper).astype(int),
        }
        for key, values in expected.items():
            if key not in frame or not np.array_equal(frame[key].to_numpy(), values.to_numpy()):
                raise ValueError("Saved input differs from audited source: " + key)
        if not audit["key_mapping_verified"]:
            raise ValueError("Source audit does not confirm key mapping")
    result = prepare_data(frame)
    before = len(result)
    if rt_min is not None:
        if not np.isfinite(rt_min) or rt_min < 0:
            raise ValueError('rt_min must be finite nonnegative seconds or None')
        result = result.loc[result.rt.ge(rt_min)].copy().reset_index(drop=True)
        validate_data(result)
    result.attrs.update(rt_min_seconds=rt_min, rows_before_rt_filter=before,
        rows_removed_by_rt_filter=before-len(result),
        rt_policy=RT_POLICY if rt_min==RT_MIN_SECONDS else 'explicit_custom_or_unfiltered')
    result.attrs.update(source_path=str(path.resolve()), source_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    return result


def smoke_subset(data, subjects=2, trials_per_cell=3):
    """Deterministic balanced small subset with existing event IDs preserved."""
    if not isinstance(subjects, int) or not isinstance(trials_per_cell, int) or subjects < 1 or trials_per_cell < 1:
        raise ValueError("Positive integer subjects and trials_per_cell required")
    data = prepare_data(data)
    ids = sorted(data.subj_idx.unique(), key=str)[:subjects]
    if len(ids) < subjects:
        raise ValueError("Not enough subjects")
    chosen = data[data.subj_idx.isin(ids)]
    groups = chosen.groupby(["subj_idx", "val", "id", "stim"], sort=True)
    if groups.size().min() < trials_per_cell:
        raise ValueError("Not enough trials in at least one design cell")
    result = groups.head(trials_per_cell).copy().reset_index(drop=True)
    validate_data(result)
    result.attrs.update(data.attrs)
    return result
