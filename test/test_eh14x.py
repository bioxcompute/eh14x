"""Basic checks that the model is reachable and behaving.

Not a benchmark. This answers "is it working", not "is it good" - those numbers
are in the README, measured on held-out taxa and held-out genes, and they sit
next to the simplest method that could produce the same answer.

What these do check is the kind of thing that breaks quietly: an endpoint that
answers with the wrong shape, a model that has silently been swapped, or a
score that no longer distinguishes real sequence from noise.

    export BIOX_API_KEY=...
    python test/test_eh14x.py
"""

from __future__ import annotations

import random
import sys

from biox import Client, BioXError

# A stretch of 18S rRNA. Real, conserved, and the sort of thing the model reads.
ECHT = ("ATGCTTGTCTCAAAGATTAAGCCATGCATGTCTCAGTATAAGCTTTTACATGGCGAAACTGCGAATGGCTC"
        "ATTAAAACAGTTATAGTTTATTTGATGGTACCTTACTACTCGGATAACCGTAGTAATTCTAGAGCTAATAC")


def geslaagd(naam: str, ok: bool, wat: str = "") -> bool:
    print(f"  {'ok  ' if ok else 'FAIL'}  {naam}" + (f"  {wat}" if wat else ""))
    return ok


def test_reachable(c: Client) -> bool:
    m = c.models()["models"][0]
    return geslaagd("the API answers and names its model", m["name"] == "EH14 X",
                    m["name"])


def test_real_beats_shuffled(c: Client) -> bool:
    """The one test worth having: real sequence should be more predictable.

    If this fails, something is wrong that no shape check would catch. The
    shuffle keeps base composition identical, so only the order differs, which
    is the only thing a sequence model can be reading.
    """
    door = list(ECHT)
    random.Random(0).shuffle(door)
    echt = c.score(ECHT, place=False, per_base=False)["nll"]
    geschud = c.score("".join(door), place=False, per_base=False)["nll"]
    return geslaagd("real sequence beats a shuffle of itself", echt < geschud,
                    f"{echt:.4f} against {geschud:.4f}")


def test_conservation_shape(c: Client) -> bool:
    r = c.score(ECHT)
    cons = r.get("conservation") or []
    return geslaagd("conservation has one value per base",
                    len(cons) == r["length"], f"{len(cons)} for {r['length']}")


def test_chimera(c: Client) -> bool:
    """The one task where this model beats the tool people actually use.

    A read, a score, and the comparison that belongs with it. If the numbers
    stop travelling with the answer, this test says so.
    """
    seq = ("ACCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGT"
           "ACGCACGGCCGGTACAGTGAAACTGCGAATGGCTCATTAAATCAGTTATCGTTTATTTGATAG"
           "TACCTTACTACATGGATAACCGTGGTAATTCTAGAGCTAATACATGCTTAAAATCTCGACCCT"
           "TTGGAAGAGATGTATTTATTAGATAAAAAATCAATGCCCTTCGGGGCTCCTTGGTGATTCATG"
           "ATAACTTCTCGAATCGCATGGCCTTGTGCTGGCGATGGTTCATTCAAATTTCTGCCCTATCA")
    r = c.chimera(seq)
    row = (r.get("results") or [{}])[0]
    ok = (row.get("call") in ("real", "chimera")
          and isinstance(r.get("auc_held_out_groups"), (int, float))
          and "9442" in str(r.get("versus_uchime", "")))
    return geslaagd("a read is called, with its AUC and the UCHIME comparison",
                    ok, f"{row.get('call')}  {row.get('score')}  "
                        f"auc {r.get('auc_held_out_groups')}")


def test_withdrawn_route_stays_withdrawn(c: Client) -> bool:
    """The variant route was withdrawn on 26 September 2026.

    EH14 X scores 0.3687 on variant effect, below the 0.5 of guessing, and two
    different models once carried one name on that route. It must answer 410
    and not quietly start working again, because a route that comes back
    without a measurement behind it is worse than one that is gone.
    """
    try:
        c._post("/v1/variant/score", {"chrom": "chr17", "pos": 43092919,
                                      "ref": "G", "alt": "A"})
    except Exception as exc:
        weg = "410" in str(exc) or "withdrawn" in str(exc).lower()
        return geslaagd("the variant route still answers 410", weg, str(exc)[:52])
    return geslaagd("the variant route still answers 410", False,
                    "it answered instead of refusing")


def test_design(c: Client) -> bool:
    doelen = c.targets()
    if not doelen:
        return geslaagd("design lists its targets", False)
    r = c.design(doelen[0]["genus"], rounds=2)
    d = r.get("designed", {})
    ok = (len(d.get("seq", "")) == 22 and set(d.get("seq", "")) <= set("ACGT")
          and r.get("tested") is False)
    return geslaagd("design returns a 22-base candidate, marked untested", ok,
                    d.get("seq", ""))


def test_design_refuses(c: Client) -> bool:
    """Designing works for a fixed list and nothing else."""
    try:
        c.design("Homo", rounds=1)
    except BioXError as exc:
        return geslaagd("design refuses a target off the list", "not one of" in str(exc))
    return geslaagd("design refuses a target off the list", False, "it answered")


if __name__ == "__main__":
    try:
        c = Client()
    except BioXError as exc:
        raise SystemExit(f"{exc}")

    tests = [test_reachable, test_real_beats_shuffled, test_conservation_shape,
             test_chimera, test_withdrawn_route_stays_withdrawn,
             test_design, test_design_refuses]
    uit = []
    for t in tests:
        try:
            uit.append(t(c))
        except Exception as exc:
            uit.append(geslaagd(t.__name__, False, f"{type(exc).__name__}: {exc}"))
    print(f"\n  {sum(uit)}/{len(uit)} passed")
    sys.exit(0 if all(uit) else 1)
