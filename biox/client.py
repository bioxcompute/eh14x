"""HTTP client for the EH14 X API. Standard library only, on purpose.

A client that needs six dependencies is a client that breaks on someone else's
machine six months from now. This uses urllib, so it works wherever Python
works.
"""

from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.request
from typing import Any

BASE = os.environ.get("BIOX_API_BASE", "https://bioxcompute.com")
AGENT = "biox-python/0.1.0"


class BioXError(RuntimeError):
    """The API said no, and the message says why."""


class Client:
    """Every call returns the parsed JSON, and carries the release that made it.

    That release string is worth keeping. Quote it in a paper and the same
    request returns the same answer later, against the same model.
    """

    def __init__(self, api_key: str | None = None, base: str = BASE,
                 timeout: float = 120.0, retries: int = 3):
        self.api_key = api_key or os.environ.get("BIOX_API_KEY")
        if not self.api_key:
            raise BioXError(
                "No API key. Request one at https://bioxcompute.com/key and set it:\n"
                "    export BIOX_API_KEY=...")
        self.base = base.rstrip("/")
        self.timeout = timeout
        self.retries = retries

    # -- the calls ---------------------------------------------------------

    def models(self) -> dict[str, Any]:
        """What the model can do, and which release you are talking to."""
        return self._get("/v1/models")

    def score(self, sequence: str, place: bool = True,
              per_base: bool = True) -> dict[str, Any]:
        """How expected is every base of a marker-gene sequence?

        The `conservation` field is the one to reach for: one number per base,
        high where the position is hard to vary. Those are the stretches that
        stayed the same across species, which is where a primer binds.
        """
        return self._post("/v1/noctua/score",
                          {"sequence": sequence, "place": place,
                           "per_base": per_base})

    def score_many(self, sequences: list[str]) -> dict[str, Any]:
        """Up to 200 sequences in one call."""
        return self._post("/v1/noctua/score", {"sequences": list(sequences)})

    def variant(self, chrom: str, pos: int, ref: str, alt: str,
                gene: str | None = None) -> dict[str, Any]:
        """How likely is this single-letter change to cause disease?

        Research use. It puts candidates in order; it does not diagnose anyone.
        """
        body = {"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt}
        if gene:
            body["gene"] = gene
        return self._post("/v1/atlas/score", body)

    def variants(self, rows: list[dict[str, Any]]) -> dict[str, Any]:
        """A whole list of variants in one call."""
        return self._post("/v1/atlas/score", {"variants": list(rows)})

    def explain(self, chrom: str, pos: int, ref: str, alt: str,
                gene: str | None = None) -> dict[str, Any]:
        """The same variant, with what the model reacted to.

        This matters more than the score. A number nobody can check is worth
        very little; this one you can argue with.
        """
        body = {"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt}
        if gene:
            body["gene"] = gene
        return self._post("/v1/atlas/explain", body)

    def generate(self, sequence: str = "", num_tokens: int = 200,
                 temperature: float = 0.9, top_k: int = 4,
                 top_p: float = 1.0, organism: str | None = None) -> dict[str, Any]:
        """Continue a sequence, or write one for a lineage."""
        body = {"sequence": sequence, "num_tokens": int(num_tokens),
                "temperature": float(temperature), "top_k": int(top_k),
                "top_p": float(top_p)}
        if organism:
            body["organism"] = organism
        return self._post("/v1/noctua/generate", body)

    def design(self, genus: str | None = None, rounds: int = 6) -> dict[str, Any]:
        """Write a detection probe for a pathogen.

        Call without `genus` to see which targets it answers for. What comes
        back is a candidate, not a validated assay, and every answer says so.
        """
        body: dict[str, Any] = {"rounds": int(rounds)}
        if genus:
            body["genus"] = genus
        return self._post("/v1/noctua/design", body)

    def targets(self) -> list[dict[str, Any]]:
        """The pathogens `design` will answer for."""
        return self.design().get("targets", [])

    # -- plumbing ----------------------------------------------------------

    def _get(self, path: str) -> dict[str, Any]:
        return self._send(urllib.request.Request(
            self.base + path, headers=self._headers()))

    def _post(self, path: str, body: dict[str, Any]) -> dict[str, Any]:
        return self._send(urllib.request.Request(
            self.base + path, json.dumps(body).encode(),
            {**self._headers(), "Content-Type": "application/json"}))

    def _headers(self) -> dict[str, str]:
        return {"Authorization": f"Bearer {self.api_key}", "User-Agent": AGENT}

    def _send(self, request) -> dict[str, Any]:
        wacht = 1.0
        for poging in range(self.retries):
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as r:
                    return json.loads(r.read())
            except urllib.error.HTTPError as exc:
                body = exc.read().decode(errors="replace")
                try:
                    melding = json.loads(body).get("error", body)
                except json.JSONDecodeError:
                    melding = body[:300]
                # 429 and 5xx are worth waiting out; a 400 will not fix itself.
                if exc.code in (429, 500, 502, 503, 504) and poging < self.retries - 1:
                    time.sleep(wacht)
                    wacht *= 2
                    continue
                raise BioXError(f"HTTP {exc.code}: {melding}") from None
            except (urllib.error.URLError, TimeoutError) as exc:
                if poging < self.retries - 1:
                    time.sleep(wacht)
                    wacht *= 2
                    continue
                raise BioXError(f"Could not reach {self.base}: {exc}") from None
        raise BioXError("unreachable")
