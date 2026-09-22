# Designed probes

Detection probes EH14 X wrote for thirteen pathogens.

| file | what |
|---|---|
| `probes.fasta` | the eight worth ordering, with their numbers in the header |
| `probes.csv` | all thirteen, including the five that did not make it and why |

Each probe is 22 bases, unmodified, orderable from any supplier as written.

## The columns

**coverage** — the share of that organism's sequences the probe binds, allowing
one mismatch. A probe with a single mismatch still binds, and a design that
ignores that is optimistic about the wrong thing.

**background_broad** — the share of *everything else* it also binds, under the
same rule, measured against 40,000 sequences from across the tree of life. Not
against the small set the probe was designed on: designing and testing on the
same data tells you nothing.

**tm_c** — melting temperature. Below 55 °C a probe tends not to bind reliably
at the temperatures a normal protocol runs at.

**orderable** — stricter than the physical filter. A strand can be physically
sound and still cover nothing, so this requires: passes the physics, covers at
least half its target, stays under one percent on the background.

## The honest part

None of these has been near a bench. Every number is computational. The only
thing that settles whether a probe works is a laboratory.

If you order one, we would like to hear what happened — and we will publish it
either way. support@bioxcompute.com
