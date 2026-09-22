"""Write a detection probe for a pathogen, and read the numbers honestly.

    python examples/03_design_a_probe.py
    python examples/03_design_a_probe.py Plasmodium
"""
import sys
from biox import Client

c = Client()

if len(sys.argv) < 2:
    print("Targets this model answers for:\n")
    for t in c.targets():
        print(f"  {t['genus']:<18} {t['disease']}")
    print("\n  python examples/03_design_a_probe.py <genus>")
    raise SystemExit(0)

r = c.design(sys.argv[1], rounds=6)
d = r["designed"]

print(f"{r['genus']} - {r['disease']}\n")
print(f"  probe        {d['seq']}  ({len(d['seq'])} bases)")
print(f"  coverage     {d['dekking']:.0%} of the target, one mismatch allowed")
print(f"  background   {d['achtergrond']:.3%} of everything else")
print(f"  Tm           {d['smeltpunt']:.1f} C")
print(f"  GC           {d['gc']:.0%}")
print(f"  usable       {d['bruikbaar']}" + ("" if d["bruikbaar"] else f" - {d['waarom']}"))

# What counting alone could find, for comparison. If the designed strand is not
# better than this, the design step did nothing and should be said so.
g = r.get("found_by_counting")
if g:
    print(f"\n  best found by counting: {g['seq']}")
    print(f"    coverage {g['dekking']:.0%}, background {g['achtergrond']:.3%}")

print(f"\n  tested in a laboratory: {r['tested']}")
print(f"  {r['use']}")
