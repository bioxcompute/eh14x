"""Place reads that match nothing in any database.

This is the part a reference database cannot do. Below five hundred metres most
of what you sequence has no named match, and a catalogue answers that with
silence. Reading answers it with the company the sequence keeps.

    python examples/04_place_unknown_reads.py reads.fasta
"""
import sys
from biox import Client


def fasta(path):
    naam, stuk = None, []
    with open(path) as f:
        for regel in f:
            regel = regel.strip()
            if regel.startswith(">"):
                if naam:
                    yield naam, "".join(stuk)
                naam, stuk = regel[1:], []
            else:
                stuk.append(regel)
    if naam:
        yield naam, "".join(stuk)


if len(sys.argv) < 2:
    raise SystemExit("usage: python examples/04_place_unknown_reads.py reads.fasta")

c = Client()
rijen = list(fasta(sys.argv[1]))
print(f"{len(rijen)} reads\n")

# Up to 200 per call; the API says so and will refuse more rather than truncate.
for i in range(0, len(rijen), 200):
    blok = rijen[i:i + 200]
    r = c.score_many([s for _, s in blok])
    for (naam, _), uit in zip(blok, r.get("results", r.get("sequences", []))):
        top = (uit.get("placement") or [{}])[0]
        groep = top.get("clazz") or top.get("phylum") or top.get("domain") or "-"
        print(f"  {naam[:34]:<36}{groep:<24}{top.get('cosine', 0):.3f}")
