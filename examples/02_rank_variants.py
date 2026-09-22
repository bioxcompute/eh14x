"""Put a list of variants of unknown significance in order.

This is the call for the queue that stalls a diagnosis: more candidates than
there are hours. It orders them so the hours go to the top.

    python examples/02_rank_variants.py
"""
from biox import Client

VARIANTS = [
    {"chrom": "chr17", "pos": 43092919, "ref": "G", "alt": "A", "gene": "BRCA1"},
    {"chrom": "chr13", "pos": 32339132, "ref": "C", "alt": "T", "gene": "BRCA2"},
    {"chrom": "chr17", "pos": 7676154, "ref": "G", "alt": "A", "gene": "TP53"},
    {"chrom": "chr7", "pos": 117559590, "ref": "G", "alt": "A", "gene": "CFTR"},
]

c = Client()
r = c.variants(VARIANTS)

print(f"{r['model']} on {r['assembly']} - {r['use']} use\n")
for v in sorted(r["variants"], key=lambda x: -x["score"]):
    print(f"  {v['score']:.4f}  {v['gene']:<8} {v['chrom']}:{v['pos']} "
          f"{v['ref']}>{v['alt']}")

# The score alone is not the useful part. Ask what the model reacted to, so a
# clinician can disagree with it on evidence rather than on faith.
top = max(r["variants"], key=lambda x: x["score"])
e = c.explain(top["chrom"], top["pos"], top["ref"], top["alt"], top.get("gene"))
print(f"\nWhat it saw for {top['gene']} {top['chrom']}:{top['pos']}:")
for a in e.get("attribution", [])[:5]:
    print(f"  offset {a['offset']:+d}  weight {a['peak_score']:.4f}")

print(f"\nMeasured: {r['measured']}")
