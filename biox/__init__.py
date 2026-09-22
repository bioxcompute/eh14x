"""A small client for the EH14 X API.

EH14 X reads biological sequence: the human genome, marker genes, and
environmental DNA from seawater and sediment. It scores how expected each base
of a sequence is, judges what a single changed letter does, places a sequence
that no database has a name for, and designs detection probes against a
pathogen.

    from biox import Client

    c = Client()                      # reads BIOX_API_KEY from the environment
    c.score("ATGCTTGTCTCAAAGATTAAGCC")
    c.variant("chr17", 43092919, "G", "A")
    c.design("Trypanosoma")

Request a key at https://bioxcompute.com/key
"""

from .client import Client, BioXError

__all__ = ["Client", "BioXError"]
__version__ = "0.1.0"
