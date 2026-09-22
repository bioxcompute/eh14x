"""Score a marker-gene sequence, and find where a primer would bind.

    python examples/01_score_a_sequence.py
"""
from biox import Client

SEQ = ("ATGCTTGTCTCAAAGATTAAGCCATGCATGTCTCAGTATAAGCTTTTACATGGCGAAACTGCGAATGGCTC"
       "ATTAAAACAGTTATAGTTTATTTGATGGTACCTTACTACTCGGATAACCGTAGTAATTCTAGAGCTAATAC")

c = Client()
r = c.score(SEQ)

print(f"{r['length']} bases, perplexity {r['perplexity']:.3f}, GC {r['gc']:.1%}")
print(f"answered by {r['model_version']}")

# Conservation is one number per base: high where the model finds that position
# hard to vary. A stretch of high conservation is a stretch that has stayed the
# same across species, which is exactly where you put a primer.
cons = r["conservation"]
venster = 22
beste = max(range(len(cons) - venster),
            key=lambda i: sum(cons[i:i + venster]) / venster)
print(f"\nMost conserved {venster}-base window starts at {beste}:")
print(f"  {SEQ[beste:beste + venster]}")
print(f"  mean conservation {sum(cons[beste:beste + venster]) / venster:.3f}")

# Placement says which reference groups it sits closest to, without looking
# anything up in a database.
for hit in r.get("placement", [])[:3]:
    naam = hit.get("clazz") or hit.get("phylum") or hit.get("domain")
    print(f"  {hit['cosine']:.3f}  {naam}")
