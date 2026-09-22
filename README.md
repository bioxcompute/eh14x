# EH14 X

A model that reads biological sequence: the human genome, marker genes, and
environmental DNA from seawater and sediment.

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

The client needs nothing beyond the standard library. `numpy` only if you run
the evaluation.

## Quickstart

```python
from biox import Client

c = Client()

# How expected is every base? conservation[] runs high where a position is
# hard to vary, which is where a primer binds.
r = c.score("ATGCTTGTCTCAAAGATTAAGCCATGCATGTCTCAGTATAAGCTTTTACATGGCGAAACT")
print(r["perplexity"], r["conservation"][:10])

# What does one changed letter do? Research use: it orders candidates.
print(c.variant("chr17", 43092919, "G", "A", gene="BRCA1"))
print(c.explain("chr17", 43092919, "G", "A"))

# Write a detection probe. What comes back is a candidate, not an assay.
print(c.design("Trypanosoma"))
```

More in [`examples/`](examples/).

## What it does

**Ranks variants nobody has ruled on.** Sequence a person and you get hundreds
of single-letter differences with no verdict attached. There are more candidates
than there are hours to spend on them, so the work is deciding which ones to
look at first. The model puts them in order and then says what it reacted to,
so a clinician can disagree on evidence.

**Places sequence that has no name.** Take a litre of seawater and sequence what
is in it. Below five hundred metres, most of what comes back matches nothing
that has been described. The model reads a sequence the way you read
handwriting, rather than looking it up, so an unnamed read still comes back with
the company it keeps.

**Writes probes that do not exist yet.** A search can only return sequences that
are already in the data. The model writes one instead: it puts the difference
where the background varies and the target does not, which is a position no
natural sequence occupies.

## Designed probes

Thirteen pathogens, every number measured, in [`probes/`](probes/) as FASTA and
CSV. Eight are worth ordering. The other five are listed with the reason they
are not, because a table where everything passes is a table nobody believes.

Each probe is 22 bases, unmodified, and can be ordered from any oligo supplier
as written.

**None has been near a bench.** Every number here is computational. The only
thing that settles whether a probe works is a laboratory. If you run one, and it
works or it does not, we would like to hear it, and we will publish the result
either way.

## How it is measured

Every number this project reports sits next to the simplest method that could
produce the same answer. A model is only worth its weights if it beats counting.

Whole groups are held out, never random rows: sequences within one gene or one
taxon are not independent of each other, and a random split lets a model answer
from memory and call it understanding. Every number carries a bootstrap
interval, because an AUC without one invites reading noise as progress.

On recognising a taxon the model has never seen, which is the question the
unnamed fraction turns on:

| | AUC | 95% interval |
|---|---|---|
| 4-mer counting | **0.774** | |
| EH14 X | 0.656 | [0.589, 0.717] |

Counting still wins, so counting is what the API serves for placement, and every
answer says which method produced it. Reproduce it with
[`evaluation/novelty.py`](evaluation/novelty.py).

## Reproducibility

Every answer carries the release that produced it:

```json
{"model": "EH14 X/1", ...}
```

Quote that string in a paper and the same request returns the same answer later,
against the same model. A release only happens when the measurement says it is
not worse than what is running. If any part of the model regresses, nothing
ships.

## Links

- [bioxcompute.com](https://bioxcompute.com), the model
- [/docs](https://bioxcompute.com/docs), API reference with curl for every call
- [/probes](https://bioxcompute.com/probes), the designed probes
- [/log](https://bioxcompute.com/log), what we tried, including what failed

## License

Apache 2.0. See [LICENSE](LICENSE).
