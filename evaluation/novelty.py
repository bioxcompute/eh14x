"""Can a model recognise a taxon it has never seen?

This is the question the dark fraction turns on. Below five hundred metres most
of what you sequence matches no named species, so a catalogue has nothing to
say. A model that reads sequence can still say "this looks like nothing I know"
without knowing what it is instead - if it can, which is what this measures.

How it is measured
------------------
Forty whole taxa are held out. Not random rows: sequences within one taxon are
not independent of each other, and a random split lets a model answer from
memory and call it recognition.

Each query is scored by mean log-likelihood per base. A known sequence should be
more predictable than an unknown one, so the measure is AUC: the chance that a
known record outranks a novel one. 0.5 is guessing.

Every model gets literally the same sequences, the same length, and no taxonomy
in the input. Comparing two models on two samples is not a comparison, and a
model that sees more context wins something that has nothing to do with the
model.

AUC and not log-likelihood, deliberately: the likelihoods of two different
models sit on different scales and must not be put next to each other. Only the
ranking is comparable.

The baseline
------------
Counting four-letter words. For each sequence, how often each of the 256 4-mers
occurs, compared against the reference groups. It is the dumbest thing that
could work, which is exactly why it belongs here: a model is only worth its
weights if it beats it.

Results, measured 22 September 2026, on the same 300 sequences of 300 bases:

    4-mer counting        0.774
    EH14 X                0.656   [0.589, 0.717]
    Evo 2 (40B)           0.540   [0.470, 0.604]

Evo 2's interval covers 0.5. Forty billion parameters do no better than chance
here, and that is not a criticism of Evo - these sequences are in no database,
so there is nothing to remember. Counting still beats both of us.

    python evaluation/novelty.py --help
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np


def auc(labels: np.ndarray, scores: np.ndarray) -> float:
    """Rank-based, so ties are handled the way they should be."""
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1)
    s = scores[order]
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + j + 2) / 2
        i = j + 1
    pos, neg = labels == 1, labels == 0
    n_p, n_n = int(pos.sum()), int(neg.sum())
    if n_p == 0 or n_n == 0:
        return float("nan")
    return (ranks[pos].sum() - n_p * (n_p + 1) / 2) / (n_p * n_n)


def interval(labels: np.ndarray, scores: np.ndarray, rounds: int = 400,
             seed: int = 0) -> tuple[float, float]:
    """Bootstrap. An AUC without one invites reading noise as progress."""
    rng = np.random.default_rng(seed)
    got = []
    for _ in range(rounds):
        pick = rng.integers(0, len(labels), len(labels))
        a = auc(labels[pick], scores[pick])
        if a == a:
            got.append(a)
    lo, hi = np.percentile(got, [2.5, 97.5])
    return float(lo), float(hi)


def score_biox(sequences: list[str]) -> list[float]:
    """Mean log-likelihood per base under EH14 X, through the public API."""
    from biox import Client
    c = Client()
    uit = []
    for s in sequences:
        r = c.score(s, place=False, per_base=False)
        # nll is the mean negative log-likelihood; flip the sign so that higher
        # means more expected, the same direction as every other scorer here.
        uit.append(-float(r["nll"]))
    return uit


def score_kmers(sequences: list[str], reference: list[str]) -> list[float]:
    """The baseline: how close is this sequence's 4-mer profile to the references?"""
    from itertools import product
    woorden = {"".join(p): i for i, p in enumerate(product("ACGT", repeat=4))}

    def profiel(s: str) -> np.ndarray:
        v = np.zeros(len(woorden))
        for i in range(len(s) - 3):
            j = woorden.get(s[i:i + 4])
            if j is not None:
                v[j] += 1
        n = v.sum()
        return v / n if n else v

    ref = np.stack([profiel(s) for s in reference])
    ref = ref / (np.linalg.norm(ref, axis=1, keepdims=True) + 1e-9)
    uit = []
    for s in sequences:
        v = profiel(s)
        v = v / (np.linalg.norm(v) + 1e-9)
        uit.append(float((ref @ v).max()))
    return uit


def run(known: list[str], novel: list[str], scorer) -> dict:
    scores = np.array(scorer(known) + scorer(novel))
    labels = np.concatenate([np.ones(len(known)), np.zeros(len(novel))])
    ok = np.isfinite(scores)
    a = auc(labels[ok], scores[ok])
    lo, hi = interval(labels[ok], scores[ok])
    return {"auc": round(float(a), 4), "ci95": [round(lo, 4), round(hi, 4)],
            "known": int(ok[:len(known)].sum()), "novel": int(ok[len(known):].sum())}


def lees_fasta(pad: str) -> list[str]:
    stukken, buf = [], []
    with open(pad) as f:
        for regel in f:
            if regel.startswith(">"):
                if buf:
                    stukken.append("".join(buf)); buf = []
            else:
                buf.append(regel.strip())
    if buf:
        stukken.append("".join(buf))
    return stukken


if __name__ == "__main__":
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    a.add_argument("--known", required=True, help="FASTA of sequences from taxa the model saw")
    a.add_argument("--novel", required=True, help="FASTA of sequences from held-out taxa")
    a.add_argument("--reference", help="FASTA of reference sequences, for the k-mer baseline")
    a.add_argument("--length", type=int, default=300, help="crop every query to this length")
    n = a.parse_args()

    known = [s[:n.length] for s in lees_fasta(n.known) if len(s) >= n.length]
    novel = [s[:n.length] for s in lees_fasta(n.novel) if len(s) >= n.length]
    print(f"{len(known)} known, {len(novel)} novel, {n.length} bases each\n")

    uit = {}
    if n.reference:
        ref = [s for s in lees_fasta(n.reference) if len(s) >= n.length]
        uit["4-mer counting"] = run(known, novel, lambda q: score_kmers(q, ref))
    if os.environ.get("BIOX_API_KEY"):
        uit["EH14 X"] = run(known, novel, score_biox)
    else:
        print("set BIOX_API_KEY to include EH14 X\n", file=sys.stderr)

    for naam, r in uit.items():
        lo, hi = r["ci95"]
        print(f"  {naam:<20}{r['auc']:.4f}   [{lo:.3f}, {hi:.3f}]")
    print(json.dumps(uit, indent=1))
