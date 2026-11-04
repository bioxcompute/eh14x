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

More in [`examples/`](examples/). A set of basic checks that the model is
reachable and behaving is in [`test/`](test/):

```bash
python test/test_eh14x.py
```

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
unnamed fraction turns on.

Two intervals are given for each number, and the second one is the one that
belongs to this task. Whole taxa are held out, so two sequences from the same
taxon are not two independent observations: resampling rows treats them as if
they were and comes out about three times too narrow. Measured on this task,
not assumed. Both are here because the narrower one is what most papers print.

| on 300 bases | AUC | over rows | over taxa |
|---|---|---|---|
| 4-mer counting | **0.776** | [0.747, 0.807] | [0.694, 0.865] |
| EH14 X | 0.753 | [0.722, 0.783] | [0.646, 0.857] |

Counting scores higher and the intervals overlap, so this measurement does not
separate the two. That is what we say now. An earlier version of this file said
"counting still wins", and its own interval did not support that: [0.722, 0.783]
contains 0.776. Over taxa the overlap is wide, and we can put a number on how
wide. The interval would have to be 4.4 times narrower before 0.776 fell outside
it.

Counting is still what the API serves for placement, because it is simpler and
nothing we have measured says to change it, and every answer says which method
produced it.

One more number that bounds the test rather than the model: 34 of the 40
held-out taxa are represented among the 400 novel sequences drawn. That is the
width of the test, and it is not the same as the width of the held-out set.

One more choice behind those numbers moves them more than the model does. The
baseline compares against one centroid per group, not against the nearest single
reference sequence: the nearest single sequence rewards a query that happens to
resemble one odd record, and gives a large group more chances to win by being
large.

The evaluation set is available on request.

## Dataset

The model reads four kinds of sequence:

| | source |
|---|---|
| human genome | GRCh38 and T2T-CHM13v2.0 |
| coding regions | GENCODE |
| reference genomes | RefSeq, one representative per species |
| environmental DNA | 18S and 16S amplicon runs from the ENA, marine and sediment |

Reference sequences come from [SILVA](https://www.arb-silva.de) (CC-BY 4.0) and
[PR2](https://pr2-database.org) (MIT), used with attribution.

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
- [/models](https://bioxcompute.com/models), every measurement, read from the reports at the moment you ask
- [/news](https://bioxcompute.com/news), the larger steps
- [probes/](probes/), the designed probes, in this repository

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
