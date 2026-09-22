# EH14 X

A model that reads the language biology is written in — the human genome, marker
genes, and environmental DNA from seawater and sediment.

It scores how expected each base of a sequence is, judges what a single changed
letter does to a person, places a sequence that no database has a name for, and
designs detection probes against a pathogen.

Access is by request: **[bioxcompute.com/request](https://bioxcompute.com/request)**

---

## Install

```bash
pip install -e .
export BIOX_API_KEY=...        # request one at bioxcompute.com/key
```

No dependencies beyond the standard library for the client. `numpy` only if you
run the evaluation.

## Quickstart

```python
from biox import Client

c = Client()

# How expected is every base? conservation[] is high where a position is hard
# to vary — which is where a primer binds.
r = c.score("ATGCTTGTCTCAAAGATTAAGCCATGCATGTCTCAGTATAAGCTTTTACATGGCGAAACT")
print(r["perplexity"], r["conservation"][:10])

# What does one changed letter do? Research use: it orders candidates.
print(c.variant("chr17", 43092919, "G", "A", gene="BRCA1"))
print(c.explain("chr17", 43092919, "G", "A"))

# Write a detection probe. What comes back is a candidate, not an assay.
print(c.design("Trypanosoma"))
```

More in [`examples/`](examples/).

## What it is for

**Variants nobody has ruled on.** Sequence a person and you get hundreds of
single-letter differences with no verdict attached. Most searches for a
diagnosis stop there — not because the work is impossible, but because there are
more candidates than hours. This puts them in order, and then says what it
reacted to, so a clinician can disagree on evidence.

**Sequence with no name.** Take a litre of seawater and sequence what is in it.
The usual method looks every read up in a catalogue, which returns nothing at
all for anything nobody has described. Below five hundred metres that is most of
what you find. This reads a sequence the way you read handwriting rather than
looking it up, so something unnamed still comes back with the company it keeps.

**Probes that do not exist yet.** Searching a database can only find sequences
that are already there. This writes one: it places the difference exactly where
the background varies and the target does not, which is a position no natural
sequence occupies.

## Designed probes

Thirteen pathogens, every number measured, in [`probes/`](probes/) as FASTA and
CSV. Eight are worth ordering; the other five are listed with the reason they
are not, because a table where everything passes is a table nobody believes.

Each is 22 bases, unmodified, and can be ordered from any oligo supplier as
written.

**None has been near a bench.** Every number here is computational, and the only
thing that settles whether a probe works is a laboratory. If you run one — and it
works, or it does not — we would like to hear it, and we will publish the result
either way.

## How it is measured

Every number this project reports sits next to the dumbest method that could
produce it. A model is only worth its weights if it beats counting.

Held-out **whole groups**, never random rows: sequences within one gene or one
taxon are not independent, and a random split lets a model answer from memory
and call it understanding. Every number carries a bootstrap interval, because an
AUC without one invites reading noise as progress.

On recognising a taxon the model has never seen — the question the dark fraction
turns on — measured on the same 300 sequences of 300 bases:

| | AUC | 95% interval |
|---|---|---|
| 4-mer counting | **0.774** | |
| EH14 X | 0.656 | [0.589, 0.717] |
| Evo 2 (40B) | 0.540 | [0.470, 0.604] |

Evo 2's interval covers 0.5: forty billion parameters do no better than chance
here. That is not a criticism of Evo — these sequences are in no database, so
there is nothing to remember. **And counting still beats both of us.** We serve
the method that wins, and we say which one it is on every answer.

Reproduce it with [`evaluation/novelty.py`](evaluation/novelty.py).

## Reproducibility

Every answer carries the release that produced it:

```json
{"model": "EH14 X/1", ...}
```

Quote that string in a paper and the same request returns the same answer later,
against the same model. A release only happens when the measurement says it is
not worse than what is running — and if any part of the model regresses, nothing
ships.

## Links

- [bioxcompute.com](https://bioxcompute.com) — the model
- [/docs](https://bioxcompute.com/docs) — API reference, with curl for every call
- [/probes](https://bioxcompute.com/probes) — the designed probes
- [/log](https://bioxcompute.com/log) — what we tried, including what failed

## License

Apache 2.0. See [LICENSE](LICENSE).
