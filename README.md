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
answer says which method produced it.

[`evaluation/novelty.py`](evaluation/novelty.py) is the code that produces that
table. It takes three files: sequences from taxa the model saw, sequences from
taxa held out of training, and a reference set with a group label per sequence.
Run it against any model you like.

**The sequences themselves are not in this repository yet.** They come from
SILVA and PR2, and we have not settled whether we may redistribute a subset, so
we would rather ship nothing than ship it wrongly. Ask and we will send the
exact set we measured on, or rebuild your own: hold out whole taxa, never random
rows, and give the baseline one centroid per group rather than the nearest
single sequence. Those two choices move the number more than the model does.

## Dataset

The model reads four kinds of sequence:

| | source |
|---|---|
| human genome | GRCh38 and T2T-CHM13v2.0 |
| coding regions | GENCODE |
| reference genomes | RefSeq, one representative per species |
| environmental DNA | 18S and 16S amplicon runs from the ENA, marine and sediment |

The environmental part is the largest by count and carries no taxonomy at all,
which is the point: a model trained only on named organisms has never seen the
thing it is meant to recognise.

Viruses that infect vertebrates are left out, along with toxin databases and
functional annotations of pathogens. What a model has never read, it cannot
write. The marker-gene corpus carries identity and not function.

Whole taxa are held out of training for the evaluation, and the list comes from
the same function that the evaluation uses to hold them out, so the two cannot
drift apart.

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

## Citation

```bibtex
@software{eh14x,
  title  = {EH14 X: a sequence model for the named and the unnamed},
  author = {{BioX}},
  year   = {2026},
  url    = {https://bioxcompute.com}
}
```

Please quote the release string that came with your answers, such as
`EH14 X/1`, so the result can be reproduced against the same model.

## License

Apache 2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE).
