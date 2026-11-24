"""Is this read a real organism, or something the PCR glued together?

This is the call where the model beats the tool people actually use. Every
answer carries its own AUC, its floor, and the comparison against UCHIME, so
the number travels with the result instead of living in a README that may be
out of date by the time you read it.

    python examples/02_chimera_or_real.py
"""
from biox import Client

# Two 18S fragments. The second one is the first 150 bases of the first, with
# 150 bases of a different lineage welded on: exactly what a PCR does when two
# half-made products find each other.
ECHT = ("ACCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGT"
        "ACGCACGGCCGGTACAGTGAAACTGCGAATGGCTCATTAAATCAGTTATCGTTTATTTGATAG"
        "TACCTTACTACATGGATAACCGTGGTAATTCTAGAGCTAATACATGCTTAAAATCTCGACCCT"
        "TTGGAAGAGATGTATTTATTAGATAAAAAATCAATGCCCTTCGGGGCTCCTTGGTGATTCATG"
        "ATAACTTCTCGAATCGCATGGCCTTGTGCTGGCGATGGTTCATTCAAATTTCTGCCCTATCA")

ANDERS = ("GTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAG"
          "CTCGTAGTTGGATTTCGGGTGGGTTCCAGCGGTCCGCCTATGGTGAGTACTGCTGTGGCCTTC"
          "CTTTCTGTCGGGGACGGGCTCCTGGGCTTCACTGTCCGGGACTCGGAGTCGGCGAGGTT")

# Both halves are cut so the welded read is exactly as long as the real one.
# If the two differed in length, length alone could be the cue, and then the
# example would be testing the wrong thing.
HALF = len(ECHT) // 2
GELAST = ECHT[:HALF] + ANDERS[:len(ECHT) - HALF]

c = Client()
r = c.chimera(sequences=[ECHT, GELAST])

print(f"{r['model']}  {r['task']}  checkpoint {r['checkpoint']}")
print(f"  AUC on held-out groups : {r['auc_held_out_groups']}")
print(f"  floor, 4-mer counting  : {r['floor_4mer_counting']}")
print(f"  against UCHIME         : {r['versus_uchime']}\n")

for naam, row in zip(("one organism", "two welded together"), r["results"]):
    if row.get("refused"):
        print(f"  {naam:<22} refused: {row['refused']}")
    else:
        print(f"  {naam:<22} {row['call']:<8} {row['score']:.4f}  "
              f"{row['bases']} bases")

# The caveat comes back with the answer, and it is the part worth reading.
# Our comparison ran at equal abundance, where both methods see only the
# letters. In a real sample UCHIME also has the abundance of each read, and we
# do not model that. So this beats UCHIME on the letters, and that is a
# narrower claim than beating UCHIME.
print(f"\n  {r['caveat']}")

# Reads outside 250 to 332 bases come back refused rather than scored. The head
# was trained on that range, and a score outside it would be a guess wearing a
# number.
kort = c.chimera("ACGT" * 20)
print(f"\n  a 80-base read: {kort['results'][0].get('refused')}")
