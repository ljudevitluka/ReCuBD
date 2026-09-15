"""Minimal Entrez client: POST requests, retries, polite rate limiting."""
import os
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL = "host-virus-refdb"


class EntrezClient:
    def __init__(self, email="", api_key_env="NCBI_API_KEY", sleep_seconds=0.4, retries=4):
        self.email = (email or "").strip()
        self.api_key = os.environ.get(api_key_env, "").strip() if api_key_env else ""
        self.sleep = float(sleep_seconds) if not self.api_key else min(float(sleep_seconds), 0.15)
        self.retries = int(retries)

    def _payload(self, params):
        p = {k: v for k, v in params.items() if v not in (None, "")}
        p["tool"] = TOOL
        if self.email:
            p["email"] = self.email
        if self.api_key:
            p["api_key"] = self.api_key
        return urllib.parse.urlencode(p).encode()

    def post(self, endpoint, **params):
        url = BASE + endpoint + ".fcgi"
        last = None
        for attempt in range(self.retries):
            try:
                req = urllib.request.Request(url, data=self._payload(params))
                with urllib.request.urlopen(req, timeout=300) as resp:
                    body = resp.read()
                time.sleep(self.sleep)
                return body
            except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
                last = exc
                time.sleep(self.sleep + 2.0 * (attempt + 1))
        raise RuntimeError("Entrez " + endpoint + " failed after " + str(self.retries) + " attempts: " + repr(last))

    def xml(self, endpoint, **params):
        return ET.fromstring(self.post(endpoint, **params))

    def esearch_count(self, db, term):
        return int(self.xml("esearch", db=db, term=term, retmax=0).findtext("Count"))

    def esearch_uids(self, db, term, page=5000):
        root = self.xml("esearch", db=db, term=term, retmax=0)
        total = int(root.findtext("Count"))
        uids = []
        for start in range(0, total, page):
            r = self.xml("esearch", db=db, term=term, retstart=start, retmax=page)
            uids.extend(e.text for e in r.findall(".//Id"))
        return uids, total

    def efetch(self, db, ids, rettype="gb", retmode="xml"):
        return self.post("efetch", db=db, id=",".join(ids), rettype=rettype, retmode=retmode)


def client_from_params(entrez_params):
    return EntrezClient(
        email=entrez_params.get("email", ""),
        api_key_env=entrez_params.get("api_key_env", "NCBI_API_KEY"),
        sleep_seconds=entrez_params.get("sleep_seconds", 0.4),
    )
