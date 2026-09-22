"""Put a list of variants of unknown significance in order.

This is the call for the queue that stalls a diagnosis: more candidates than
there are hours. It orders them so the hours go to the top.

    python examples/02_rank_variants.py
"""
from biox import Client

VARIANTS = [
    {"chrom": "chr17", "pos": 43092919, "ref": "G", "alt": "A", "gene": "BRCA1"},
    {"chrom": "chr13", "pos": 32339132, "ref": "G", "alt": "A", "gene": "BRCA2"},
    {"chrom": "chr17", "pos": 7676154, "ref": "G", "alt": "A", "gene": "TP53"},
    {"chrom": "chr7", "pos": 117559590, "ref": "A", "alt": "G", "gene": "CFTR"},
    # Deliberately wrong, to show what happens. The reference base at this
    # position is not a T, and the API says so instead of scoring it anyway.
    {"chrom": "chr17", "pos": 43092919, "ref": "T", "alt": "C", "gene": "BRCA1"},
]

c = Client()
r = c.variants(VARIANTS)

print(f"{r['model']} on {r['assembly']} - {r['use']} use\n")

gescoord = [v for v in r["variants"] if "score" in v]
mislukt = [v for v in r["variants"] if "error" in v]

for v in sorted(gescoord, key=lambda x: -x["score"]):
    print(f"  {v['score']:.4f}  {v['gene']:<8} {v['chrom']}:{v['pos']} "
          f"{v['ref']}>{v['alt']}")

# A variant whose reference base does not match GRCh38 is refused rather than
# scored. That check is worth more than it looks: it catches a wrong assembly
# or an off-by-one coordinate before it turns into a result you believe.
if mislukt:
    print(f"\n  {len(mislukt)} refused:")
    for v in mislukt:
        print(f"    {v['gene']} {v['chrom']}:{v['pos']} - {v['error']}")

# The score alone is not the useful part. Ask what the model reacted to, so a
# clinician can disagree with it on evidence rather than on faith.
if gescoord:
    top = max(gescoord, key=lambda x: x["score"])
    e = c.explain(top["chrom"], top["pos"], top["ref"], top["alt"], top.get("gene"))
    print(f"\nWhat it saw for {top['gene']} {top['chrom']}:{top['pos']}:")
    for a in e.get("attribution", [])[:5]:
        print(f"  offset {a['offset']:+d}  weight {a['peak_score']:.4f}")

print(f"\nMeasured: {r['measured']}")
