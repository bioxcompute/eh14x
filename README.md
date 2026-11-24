# EH14 X

A model that reads biological sequence: the human genome, marker genes, and
environmental DNA from seawater and sediment.

It scores how expected each base of a sequence is, tells a real read from a
PCR artefact, places a sequence that no database has a name for, and designs
detection probes against a pathogen.

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

# Real organism, or two things the PCR welded together? The head was trained
# on 250 to 332 bases and refuses a read outside that range. Every answer
# carries its own AUC, its floor, and the comparison against UCHIME.
q = c.chimera(read)                     # one 18S read, 250 to 332 bases
print(q["results"][0]["call"], q["versus_uchime"])

# Write a detection probe. What comes back is a candidate, not an assay.
print(c.design("Trypanosoma"))
```

The variant route this client used to carry was withdrawn on 26 September
2026. EH14 X scores 0.3687 on variant effect, below the 0.5 of guessing, and
two different models once answered on that one name. It returns 410 rather than
404, because it existed and was taken away on purpose.

More in [`examples/`](examples/). A set of basic checks that the model is
reachable and behaving is in [`test/`](test/):

```bash
python test/test_eh14x.py
```

## What it does

**Tells a real read from a PCR artefact.** Amplify DNA out of seawater and some
of what comes back is two half-made products welded together. Those welds look
like new organisms, and a pipeline that keeps them reports diversity that is not
there. This is the one task here where the model beats the tool people actually
use: 0.9735 against 0.9442 for UCHIME in vsearch 2.32 on the same set, a
difference of 0.0293 with interval [0.0208, 0.0385].

That comparison ran at equal abundance and abskew 1.0, so both methods saw only
the letters. In a real sample UCHIME also has the abundance of every read, and
we do not model it, so this beats UCHIME on the letters and that is a narrower
claim than beating UCHIME. Our own 4-mer baseline scored 0.4989 here, which
means that baseline was too easy and the honest comparison is the one above.

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

| on 512 bases | AUC | over rows | over taxa |
|---|---|---|---|
| 4-mer counting | **0.823** | [0.796, 0.853] | [0.746, 0.901] |
| EH14 X | 0.797 | [0.765, 0.825] | [0.692, 0.894] |

More context is easier for both methods, so a model read at one length may only
be compared with a baseline measured at that same length. Both lengths are here
because the published comparisons against other models were run at 300 bases.

Counting scores higher at both lengths and the intervals overlap at both, so
this measurement does not separate the two. That is what we say now. An earlier
version of this file said "counting still wins", and its own interval did not
support that: [0.722, 0.783] contains 0.776.

How wide is the overlap? Over taxa the interval would have to be 4.4 times
narrower at 300 bases, and 3.8 times narrower at 512, before the counting score
fell outside it. Over rows it would take only 1.3 and 1.1, which is the other
reason the unit matters: on the wrong unit this looks like a near miss.

A sharper test than two separate intervals is possible here, and it was not
being run. Both methods see exactly the same sequences in the same order, which
is checked per position rather than assumed, so the difference itself can be
resampled. Two separate intervals each carry the full sampling noise, and that
noise is shared between the two methods; resampling the difference removes it.
Paired over taxa, the difference is -0.024 [-0.076, +0.026] at 300 bases and
-0.026 [-0.092, +0.035] at 512, and the model scores higher in about a quarter
of 2000 resamples. Zero sits inside both intervals, so the sharper test says the
same thing as the wider one: this measurement does not separate the two.

The floor of that test was measured too, because a test that always finds
something finds nothing. Two genuinely different baselines whose scores differ
by 0.001 give a difference interval that contains zero, and one that differs by
0.005 does not, so the resolution on this task is about 0.005. All of these
numbers are at 400 sequences per side for both methods: a baseline measured on a
different number of sequences has a different interval width, and the width is
what decides whether anything is separated.

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

The evaluation set is published. **[dark-taxa-novelty
v1.0](https://github.com/bioxcompute/dark-taxa-novelty)** is the frozen split:
1,760 sequences, forty whole taxa held out, CC-BY 4.0, with the baseline in one
numpy file so the row above can be reproduced before anyone tries to beat it.
Measured there, the baseline reaches 0.7741 on full reads, 0.7578 at 512 bases
and 0.7329 at 300. Those are its own forty held-out taxa and not the forty in
the table above, so the two sets of numbers are not interchangeable.

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
  organization = {bioX},
  author = {Snooij, Niels Karel},
  year   = {2026},
  url    = {https://bioxcompute.com}
}
```

Please quote the release string that came with your answers, such as
`EH14 X/1`, so the result can be reproduced against the same model.

## License

Apache 2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE).
