#!/usr/bin/env python3
# SPDX-License-Identifier: CC BY-NC 4.0
# Copyright (c) 2026 Daniel Pérez
"""
sonify_tissuespectf.py

Create one MIDI composition per TissueSpectF state.

Concept
-------
- core invariant spectral peaks = shared accompaniment
- condition consensus peaks     = state-specific melody

Biological-to-musical mapping
-----------------------------
Spectral frequency f=k/N:
    mapped logarithmically to MIDI pitch.
    Higher genomic spectral frequency -> higher musical pitch.

Phase:
    mapped to a small sub-beat timing displacement.

Peak strength:
    invariant_score / meta_score / final_power -> MIDI velocity.

Period:
    mapped inversely to note duration category:
    long genomic periods -> longer notes.

Evidence class:
    robust    -> stronger melodic articulation
    candidate -> softer articulation

This is a sonification: musical choices are deterministic and mapping tables are
written alongside every MIDI so the transformation remains inspectable.

Requirements
------------
pip install pandas numpy mido

Example
-------
python3 scripts/sonify_tissuespectf.py \
  --library-dir results/library_domains \
  --conditions Normal_histology,F0,F1,F2,F3,F4

--library-dir defaults to $TSF_LIBRARY_DIR, then to
$TSF_RESULTS_DIR/library_domains, then to <repo>/results/library_domains, the
same resolution the rest of the pipeline uses. --out-dir defaults to a
"sonification" folder inside it.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    from mido import Message, MetaMessage, MidiFile, MidiTrack, bpm2tempo
except ImportError as exc:
    raise SystemExit(
        "Missing dependency 'mido'. Install with:\n"
        "  pip install mido pandas numpy"
    ) from exc


VERSION = "2026-09-01-sonification-v1"

TICKS_PER_BEAT = 480
BEATS_PER_BAR = 4
DEFAULT_BARS = 16
DEFAULT_BPM = 82

# A restrained modal palette. Frequencies still determine pitch height;
# snapping only makes the result listenable.
# D Dorian pitch classes: D E F G A B C
D_DORIAN = {0, 2, 4, 5, 7, 9, 10}  # relative to C pitch classes


def default_library_dir() -> str:
    """Resolve the library the same way the R pipeline resolves its paths.

    The default used to be "results_gencode_v3/library_domains", a directory
    that existed on one machine. A default path is a statement about layout,
    never about a location.
    """
    explicit = os.environ.get("TSF_LIBRARY_DIR", "")
    if explicit:
        return explicit
    results = os.environ.get("TSF_RESULTS_DIR", "")
    if not results:
        results = os.path.join(os.environ.get("TSF_ROOT", os.getcwd()), "results")
    return os.path.join(results, "library_domains")


def parse_args():
    p = argparse.ArgumentParser(
        description="Sonify TissueSpectF condition spectra into MIDI.")
    p.add_argument(
        "--library-dir",
        default=default_library_dir(),
        help="Condition library written by build_final_condition_spectra.R "
             "(default: $TSF_LIBRARY_DIR, else $TSF_RESULTS_DIR/library_domains).",
    )
    p.add_argument(
        "--out-dir",
        default=None,
    )
    p.add_argument(
        "--conditions",
        default="Normal_histology,F0,F1,F2,F3,F4",
    )
    p.add_argument("--bpm", type=int, default=DEFAULT_BPM)
    p.add_argument("--bars", type=int, default=DEFAULT_BARS)
    p.add_argument(
        "--max-invariants",
        type=int,
        default=16,
        help="Top core invariant modes used in accompaniment.",
    )
    p.add_argument(
        "--max-condition-peaks",
        type=int,
        default=16,
        help="Maximum consensus peaks used in each melody.",
    )
    p.add_argument(
        "--candidate-velocity-scale",
        type=float,
        default=0.72,
    )
    return p.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def first_existing(df: pd.DataFrame, cols: Iterable[str], default=None):
    for c in cols:
        if c in df.columns:
            return c
    return default


def norm01(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    finite = np.isfinite(x)
    out = pd.Series(np.zeros(len(x), dtype=float), index=x.index)
    if not finite.any():
        return out

    lo = float(x[finite].min())
    hi = float(x[finite].max())
    if math.isclose(lo, hi):
        out.loc[finite] = 0.65
    else:
        out.loc[finite] = (x.loc[finite] - lo) / (hi - lo)
    return out.clip(0, 1)


def spectral_frequency(df: pd.DataFrame) -> pd.Series:
    if "freq" in df.columns:
        f = pd.to_numeric(df["freq"], errors="coerce")
    elif {"k", "N"}.issubset(df.columns):
        f = (
            pd.to_numeric(df["k"], errors="coerce")
            / pd.to_numeric(df["N"], errors="coerce")
        )
    else:
        raise ValueError("Peak table needs freq or both k and N.")
    return f


def spectral_period(df: pd.DataFrame) -> pd.Series:
    if "period" in df.columns:
        return pd.to_numeric(df["period"], errors="coerce")
    f = spectral_frequency(df)
    return 1.0 / f.replace(0, np.nan)


def snap_to_d_dorian(note: int, low: int, high: int) -> int:
    note = max(low, min(high, int(round(note))))
    if note % 12 in D_DORIAN:
        return note

    for d in range(1, 7):
        up = note + d
        down = note - d
        if up <= high and up % 12 in D_DORIAN:
            return up
        if down >= low and down % 12 in D_DORIAN:
            return down
    return note


def freq_to_pitch(
    f: float,
    fmin: float,
    fmax: float,
    low: int,
    high: int,
) -> int:
    f = max(float(f), 1e-12)
    fmin = max(float(fmin), 1e-12)
    fmax = max(float(fmax), fmin * 1.000001)

    z = (math.log2(f) - math.log2(fmin)) / (
        math.log2(fmax) - math.log2(fmin)
    )
    raw = low + z * (high - low)
    return snap_to_d_dorian(int(round(raw)), low, high)


def phase_to_fraction(phase: float) -> float:
    if not np.isfinite(phase):
        return 0.0
    return ((float(phase) + math.pi) % (2 * math.pi)) / (2 * math.pi)


def duration_from_period(period: float, pmin: float, pmax: float) -> float:
    if not np.isfinite(period):
        return 0.5
    if math.isclose(pmin, pmax):
        return 1.0

    z = (math.log2(period) - math.log2(pmin)) / (
        math.log2(pmax) - math.log2(pmin)
    )
    # 0.25 to 1.5 beats
    return float(np.clip(0.25 + 1.25 * z, 0.25, 1.5))


def to_tick(beats: float) -> int:
    return max(0, int(round(beats * TICKS_PER_BEAT)))


def add_note_events(
    track: MidiTrack,
    events: list[dict],
    channel: int,
):
    midi_events = []
    for ev in events:
        start = to_tick(ev["start_beat"])
        end = to_tick(ev["start_beat"] + ev["duration_beats"])
        if end <= start:
            end = start + max(1, TICKS_PER_BEAT // 8)

        midi_events.append(
            (start, 1, Message(
                "note_on",
                note=int(ev["midi_note"]),
                velocity=int(ev["velocity"]),
                channel=channel,
                time=0,
            ))
        )
        midi_events.append(
            (end, 0, Message(
                "note_off",
                note=int(ev["midi_note"]),
                velocity=0,
                channel=channel,
                time=0,
            ))
        )

    # At the same tick, note_off before note_on.
    midi_events.sort(key=lambda x: (x[0], x[1]))

    prev = 0
    for tick, _, msg in midi_events:
        msg.time = tick - prev
        track.append(msg)
        prev = tick


def prepare_invariants(path: Path, max_peaks: int) -> pd.DataFrame:
    df = read_tsv(path).copy()

    if "invariant_class" in df.columns:
        core = df[df["invariant_class"].astype(str) == "core_invariant"].copy()
        if len(core):
            df = core

    score_col = first_existing(
        df,
        ["invariant_score", "median_power", "power"],
    )
    if score_col is None:
        df["_strength"] = 1.0
    else:
        df["_strength"] = pd.to_numeric(df[score_col], errors="coerce")

    phase_col = first_existing(
        df,
        ["mean_phase", "mean_phase_between_cohorts", "phase"],
    )
    df["_phase"] = (
        pd.to_numeric(df[phase_col], errors="coerce")
        if phase_col
        else 0.0
    )
    df["_freq"] = spectral_frequency(df)
    df["_period"] = spectral_period(df)

    df = df[np.isfinite(df["_freq"]) & (df["_freq"] > 0)].copy()
    df = df.sort_values("_strength", ascending=False).head(max_peaks)
    return df.reset_index(drop=True)


def prepare_condition(path: Path, max_peaks: int) -> pd.DataFrame:
    df = read_tsv(path).copy()

    # Only selected signature rows should normally be present, but protect
    # against background rows if a broader table is supplied.
    if "signature_class" in df.columns:
        allowed = df["signature_class"].astype(str).isin(["robust", "candidate"])
        if allowed.any():
            df = df[allowed].copy()

    score_col = first_existing(
        df,
        ["meta_score", "condition_specific_power", "final_power", "power"],
    )
    if score_col is None:
        df["_strength"] = 1.0
    else:
        df["_strength"] = pd.to_numeric(df[score_col], errors="coerce")

    phase_col = first_existing(
        df,
        ["mean_phase_between_cohorts", "mean_phase", "phase"],
    )
    df["_phase"] = (
        pd.to_numeric(df[phase_col], errors="coerce")
        if phase_col
        else 0.0
    )
    df["_freq"] = spectral_frequency(df)
    df["_period"] = spectral_period(df)

    if "signature_class" not in df.columns:
        df["signature_class"] = "selected"

    df = df[np.isfinite(df["_freq"]) & (df["_freq"] > 0)].copy()

    # Evidence first, then strength.
    class_rank = {"robust": 0, "candidate": 1, "selected": 2}
    df["_class_rank"] = (
        df["signature_class"]
        .astype(str)
        .map(class_rank)
        .fillna(3)
    )

    df = df.sort_values(
        ["_class_rank", "_strength"],
        ascending=[True, False],
    ).head(max_peaks)

    return df.reset_index(drop=True)


def accompaniment_events(
    inv: pd.DataFrame,
    total_beats: int,
    global_fmin: float,
    global_fmax: float,
) -> tuple[list[dict], pd.DataFrame]:
    if inv.empty:
        return [], pd.DataFrame()

    strength = norm01(inv["_strength"])
    pmin = max(float(inv["_period"].min()), 1e-9)
    pmax = max(float(inv["_period"].max()), pmin)

    # The accompaniment is a repeating arpeggio. Each invariant peak appears
    # once per cycle; biological phase offsets it within its slot.
    cycle_beats = 8.0
    slots = max(1, len(inv))
    slot_width = cycle_beats / slots

    events = []
    rows = []

    for i, row in inv.iterrows():
        pitch = freq_to_pitch(
            row["_freq"], global_fmin, global_fmax, low=36, high=62
        )
        vel = int(round(32 + 40 * strength.iloc[i]))
        dur = duration_from_period(row["_period"], pmin, pmax)
        phase_fraction = phase_to_fraction(row["_phase"])

        base_offset = i * slot_width
        phase_offset = phase_fraction * min(slot_width * 0.45, 0.35)

        start = base_offset + phase_offset
        while start < total_beats:
            ev = {
                "start_beat": start,
                "duration_beats": min(dur, max(0.25, slot_width * 1.5)),
                "midi_note": pitch,
                "velocity": vel,
            }
            events.append(ev)
            rows.append({
                "role": "invariant_accompaniment",
                "chr": row.get("chr", ""),
                "N": row.get("N", ""),
                "k": row.get("k", ""),
                "freq": row["_freq"],
                "period": row["_period"],
                "phase": row["_phase"],
                "source_strength": row["_strength"],
                "midi_note": pitch,
                "velocity": vel,
                "start_beat": start,
                "duration_beats": ev["duration_beats"],
                "invariant_class": row.get("invariant_class", "core_invariant"),
            })
            start += cycle_beats

    return events, pd.DataFrame(rows)


def melody_events(
    cond: pd.DataFrame,
    total_beats: int,
    global_fmin: float,
    global_fmax: float,
    candidate_scale: float,
) -> tuple[list[dict], pd.DataFrame]:
    if cond.empty:
        return [], pd.DataFrame()

    strength = norm01(cond["_strength"])
    pmin = max(float(cond["_period"].min()), 1e-9)
    pmax = max(float(cond["_period"].max()), pmin)

    # Distribute selected peaks through the full composition in spectral order.
    # Phase gives a small within-slot displacement.
    melody = cond.sort_values("_freq").reset_index(drop=True)
    strength_lookup = norm01(melody["_strength"])

    slot_width = total_beats / max(1, len(melody))
    events = []
    rows = []

    for i, row in melody.iterrows():
        pitch = freq_to_pitch(
            row["_freq"], global_fmin, global_fmax, low=60, high=88
        )

        cls = str(row.get("signature_class", "selected"))
        evidence_scale = candidate_scale if cls == "candidate" else 1.0

        vel = int(round((60 + 55 * strength_lookup.iloc[i]) * evidence_scale))
        vel = int(np.clip(vel, 30, 120))

        dur = duration_from_period(row["_period"], pmin, pmax)
        phase_fraction = phase_to_fraction(row["_phase"])
        phase_offset = (phase_fraction - 0.5) * min(slot_width * 0.40, 0.45)

        start = i * slot_width + max(-0.35, phase_offset)
        start = max(0.0, min(total_beats - 0.1, start))

        ev = {
            "start_beat": start,
            "duration_beats": min(dur, max(0.25, slot_width * 0.8)),
            "midi_note": pitch,
            "velocity": vel,
        }
        events.append(ev)

        rows.append({
            "role": "condition_melody",
            "chr": row.get("chr", ""),
            "N": row.get("N", ""),
            "k": row.get("k", ""),
            "freq": row["_freq"],
            "period": row["_period"],
            "phase": row["_phase"],
            "source_strength": row["_strength"],
            "signature_class": cls,
            "midi_note": pitch,
            "velocity": vel,
            "start_beat": start,
            "duration_beats": ev["duration_beats"],
        })

    return events, pd.DataFrame(rows)


def create_midi(
    condition: str,
    inv_events: list[dict],
    melody: list[dict],
    path: Path,
    bpm: int,
):
    mid = MidiFile(ticks_per_beat=TICKS_PER_BEAT)

    meta = MidiTrack()
    mid.tracks.append(meta)
    meta.append(MetaMessage(
        "track_name",
        name=f"TissueSpectF {condition}",
        time=0,
    ))
    meta.append(MetaMessage(
        "set_tempo",
        tempo=bpm2tempo(bpm),
        time=0,
    ))
    meta.append(MetaMessage(
        "time_signature",
        numerator=4,
        denominator=4,
        time=0,
    ))

    accomp = MidiTrack()
    mid.tracks.append(accomp)
    accomp.append(MetaMessage(
        "track_name",
        name="Core invariant accompaniment",
        time=0,
    ))
    # Acoustic Grand Piano, GM program 0
    accomp.append(Message("program_change", program=0, channel=0, time=0))
    add_note_events(accomp, inv_events, channel=0)

    mel = MidiTrack()
    mid.tracks.append(mel)
    mel.append(MetaMessage(
        "track_name",
        name=f"{condition} consensus melody",
        time=0,
    ))
    # Vibraphone, GM program 11
    mel.append(Message("program_change", program=11, channel=1, time=0))
    add_note_events(mel, melody, channel=1)

    path.parent.mkdir(parents=True, exist_ok=True)
    mid.save(path)


def main():
    args = parse_args()

    library_dir = Path(args.library_dir)
    out_dir = (
        Path(args.out_dir)
        if args.out_dir
        else library_dir / "sonification"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    conditions = [
        x.strip()
        for x in args.conditions.split(",")
        if x.strip()
    ]

    invariant_path = library_dir / "invariants" / "invariant_core.tsv"
    if not invariant_path.exists():
        fallback = library_dir / "invariants" / "invariant_spectrum.tsv"
        if fallback.exists():
            invariant_path = fallback
        else:
            raise SystemExit(
                f"Could not find invariant_core.tsv or invariant_spectrum.tsv "
                f"under {library_dir / 'invariants'}"
            )

    inv = prepare_invariants(invariant_path, args.max_invariants)

    condition_tables = {}
    for condition in conditions:
        path = library_dir / f"condition_signature_{condition}.tsv"
        if not path.exists():
            print(f"WARNING: missing {path}; skipping {condition}")
            continue
        condition_tables[condition] = prepare_condition(
            path,
            args.max_condition_peaks,
        )

    if not condition_tables:
        raise SystemExit("No condition signature tables found.")

    all_freq = [inv["_freq"]]
    all_freq.extend(df["_freq"] for df in condition_tables.values() if len(df))
    all_freq = pd.concat(all_freq, ignore_index=True)
    all_freq = all_freq[np.isfinite(all_freq) & (all_freq > 0)]

    global_fmin = float(all_freq.min())
    global_fmax = float(all_freq.max())

    total_beats = args.bars * BEATS_PER_BAR

    inv_events, inv_map = accompaniment_events(
        inv,
        total_beats,
        global_fmin,
        global_fmax,
    )

    inv_map.to_csv(
        out_dir / "invariant_accompaniment_mapping.tsv",
        sep="\t",
        index=False,
    )

    manifest_rows = []

    for condition, cond in condition_tables.items():
        melody, melody_map = melody_events(
            cond,
            total_beats,
            global_fmin,
            global_fmax,
            args.candidate_velocity_scale,
        )

        midi_path = out_dir / f"TissueSpectF_{condition}.mid"
        create_midi(
            condition,
            inv_events,
            melody,
            midi_path,
            args.bpm,
        )

        melody_map.to_csv(
            out_dir / f"TissueSpectF_{condition}_melody_mapping.tsv",
            sep="\t",
            index=False,
        )

        manifest_rows.append({
            "condition": condition,
            "midi_file": midi_path.name,
            "n_invariant_modes": len(inv),
            "n_condition_peaks": len(cond),
            "n_robust": int(
                (cond["signature_class"].astype(str) == "robust").sum()
            ) if "signature_class" in cond.columns else 0,
            "n_candidate": int(
                (cond["signature_class"].astype(str) == "candidate").sum()
            ) if "signature_class" in cond.columns else 0,
            "bpm": args.bpm,
            "bars": args.bars,
            "version": VERSION,
        })

        print(
            f"{condition:20s} -> {midi_path.name} "
            f"(melody peaks={len(cond)})"
        )

    pd.DataFrame(manifest_rows).to_csv(
        out_dir / "sonification_manifest.tsv",
        sep="\t",
        index=False,
    )

    readme = f"""TissueSpectF spectral sonification
Version: {VERSION}

Shared accompaniment
--------------------
Source: {invariant_path}
Top invariant modes used: {len(inv)}
Role: same accompaniment in every state.

Condition melody
----------------
Source: condition_signature_<condition>.tsv
Robust peaks have full articulation.
Candidate peaks are attenuated by factor {args.candidate_velocity_scale}.

Mapping
-------
spectral frequency (k/N) -> logarithmic MIDI pitch
phase                    -> within-slot timing displacement
peak strength            -> velocity
period                    -> duration
evidence class            -> articulation strength

Musical constraints
-------------------
Pitches are snapped to D Dorian only after their spectral height is calculated.
This improves listenability while preserving the rank/order of genomic
spectral frequency.

Important interpretation
------------------------
The MIDI is a deterministic sonification of the spectral representation.
It is not evidence that biological systems literally encode musical notes.
The TSV mapping files are the auditable bridge between every spectral peak
and every MIDI event.
"""
    (out_dir / "README_sonification.txt").write_text(readme)

    print(f"\nWrote sonifications to: {out_dir}")


if __name__ == "__main__":
    main()
