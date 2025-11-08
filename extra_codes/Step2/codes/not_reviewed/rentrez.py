#!/usr/bin/env python3
"""
Classify genome assemblies into host / source categories
--------------------------------------------------------

• Fetches BioSample metadata one-by-one with a 0.5-s delay (retries up to 5× on HTTP 429).
• Robust XML parsing (Bio.Entrez → ElementTree fallback).
• Host clues are searched in BOTH “host” and “isolation_source”.
• Food words (shrimp, chicken, turkey …) trigger “Food” only when they come
  from isolation_source, not from host.
• Whole-word matching avoids pig → pigeon mix-ups.
• Fallback labels: “Gut (Undefined host)”, “Human (Undefined tissue)”, etc.
"""

import argparse, io, os, re, time, xml.etree.ElementTree as ET
from urllib.error import HTTPError

import pandas as pd
from Bio import Entrez


# ───────────────────────── CLI ──────────────────────────
cli = argparse.ArgumentParser(description="Classify genome accessions")
cli.add_argument("--input_csv", required=True)
cli.add_argument("--accession_column", default="Genome")
cli.add_argument("--output_csv", default="classified_genomes.csv")
cli.add_argument("--api_key", help="NCBI API key or env NCBI_API_KEY")
args = cli.parse_args()

# ───────────── Entrez settings & throttle ───────────────
Entrez.email = "felipe.morillo@crg.eu"
Entrez.api_key = args.api_key or os.getenv("NCBI_API_KEY")
RATE, MAX_TRY = 0.5, 5           # 0.5 s between calls, 5 retries on 429


def _efetch(db: str, _id: str, **kw):
    """EFetch with delay + exponential back-off for HTTP 429."""
    for attempt in range(MAX_TRY):
        try:
            time.sleep(RATE)
            return Entrez.efetch(db=db, id=_id, **kw)
        except HTTPError as err:
            if err.code == 429:
                time.sleep((2 ** attempt) * RATE)
            else:
                raise
    raise RuntimeError(f"Too many HTTP 429 errors for {_id}")


# ───────────────── XML → dict helper ────────────────────
def _parse_bs(handle):
    xml = handle.read(); handle.close()
    try:
        rec = Entrez.read(io.StringIO(xml))
        attrs = rec["BioSampleSet"]["BioSample"][0]["Attributes"]
        return {a["attribute_name"].lower(): a["content"]
                for a in attrs if "attribute_name" in a}
    except Exception:                       # ElementTree fallback
        root, meta = ET.fromstring(xml), {}
        for node in root.iter():
            if node.tag.endswith("Attribute"):
                name = (node.attrib.get("attribute_name")
                        or node.attrib.get("harmonized_name")
                        or node.attrib.get("display_name"))
                if name:
                    meta[name.lower()] = (node.text or "")
        return meta


# ──────────────── vocabulary sets ────────────────
GUT   = {"gut","gastrointestinal","digestive","intestine","intestinal","cecum",
         "caecum","ceacal","caecal","cecal","hindgut","ileum","jejunum","duodenum",
         "colon","rectum","rectal","stool","feces","faeces","fecal","rumen", "lumen", "luminal",
         "ruminal","gizzard","cecal content","coprolite","midgut", "colonic", "mucusa"}
OTHER_ORG   = {"hepatopancreas"}         
ORAL  = {"oral","mouth","saliva","dental","plaque","supragingival",
         "subgingival","tongue","buccal","Salivary_gland","salivary"}
RESP  = {"nasal","nares","nose","pharyngeal","throat","sputum",
         "lung","bronchus","airway","respiratory","trachea"}
URO   = {"vagina","vaginal","cervix","urogenital","urine","bladder",
         "uterus","eggs","egg","semen"}
SKIN  = {"skin","epidermis","dermis","cutaneous","sebaceous","pore"}
ENV   = {"soil","mud","water","marine","basin","freshwater","river","lake",
         "sludge","pond","aquaculture","sediment","environment","swab",
         "mat","forest","petroleum"}
FOOD  = {"milk","cheese","yogurt","ferment","food","salami","sausage","shrimp",
         "turkey","chicken","kimchi","kombucha","beef","fiber",
         "comminuted","feed","silage","wine"}
PLANT = {"plant","root","leaf","stem","seed","rhizosphere","phyllosphere",
         "endophyte"}

MOUSE   = {"mus","mouse","mus musculus","murine","murinae","muridae","mice"}
RAT     = {"rat","rattus","rattus norvegicus","rats"}
RODENTS = {"hamster","apodemus","mesocricetus","cricetulus","gerbil","vole",
           "squirrel","guinea pig","cavia","chinchilla","rodent","myodes",
           "marmot","wood mouse","yellow-necked mouse", "apodemus","castor"}
HUMAN   = {"human","homo","homo sapiens","patient","h. sapiens"}
MAMMAL  = {"pig","sus","boar","cow","bos","bovine","cattle","ruminant","sheep","ovis","goat",
           "capra","deer","cervus","yak","canis","dog","cat","felis","rabbit",
           "oryctolagus","camel","horse","equus","phascolarctos","koala",
           "wombat","vombatus","bison","water buffalo","macaca","baboon","Myotragus"}
OTHER   = {"chicken","gallus","duck","anas","fish","danio","zebrafish","trout",
           "salmo","tilapia","oreochromis","frog","xenopus","amphibian",
           "reptile","snake","serpentes","lizard","shrimp","bird","aves",
           "pigeon","pigeons","papio","penaeus","tachysurus","gasterosteus"}


# ──────────── whole-word matcher ────────────
def _any(keys, text: str) -> bool:
    """Return True if *any* keyword appears as a whole word (case-insensitive)."""
    t = text.lower()
    for k in keys:
        if re.search(rf"\b{re.escape(k.lower())}\b", t):
            return True
    return False


# ─────── site detector (isolation_source only) ───────
def site_category(iso_text: str) -> str | None:
    t = iso_text.lower()
    if _any(GUT,   t):   return "gut"
    if _any(ORAL,  t) or _any(RESP, t) or _any(URO, t) or _any(SKIN, t) or _any(OTHER_ORG, t):
        return "non-gut"
    if _any(ENV,   t):   return "environment"
    if _any(FOOD,  t):   return "food"
    if _any(PLANT, t):   return "plant"
    return None


# ─────────────── classifier ───────────────
def classify_source(host_field: str, iso_field: str) -> str:
    """Return final 'Type' label for one sample."""
    host_group, host_origin = None, None
    for name, keyset in [
        ("Rat", RAT), ("Mouse", MOUSE), ("Other rodents", RODENTS),
        ("Human", HUMAN), ("Other mammals", MAMMAL), ("Other animals", OTHER)
    ]:
        if _any(keyset, host_field):
            host_group, host_origin = name, "host"
            break
        if _any(keyset, iso_field):
            host_group, host_origin = name, "iso"
            break

    site = site_category(iso_field)

    if site == "gut":
        return f"{host_group} gut" if host_group else "Gut (Undefined host)"

    if site == "non-gut":
        return f"{host_group} non-gut" if host_group else "Non-gut (Undefined host)"

    if site == "food" and host_origin != "host":
        return "Food"

    if site in {"environment", "plant"}:
        return site.capitalize()

    if host_group:
        return f"{host_group} (Undefined tissue)"

    return "Undefined"


# ── fetch ▪︎ parse ▪︎ classify a single accession ──
def classify_accession(acc: str) -> dict:
    try:
        # 1 ▸ assembly UID
        time.sleep(RATE)
        uid_list = Entrez.read(Entrez.esearch(db="assembly", term=acc))["IdList"]
        if not uid_list:
            return {"Genome": acc, "Isolation source": "NA", "Type": "Undefined"}
        uid = uid_list[0]

        # 2 ▸ BioSample accession
        time.sleep(RATE)
        bs_acc = Entrez.read(Entrez.esummary(db="assembly", id=uid)) \
                       ["DocumentSummarySet"]["DocumentSummary"][0] \
                       .get("BioSampleAccn")
        if not bs_acc:
            return {"Genome": acc, "Isolation source": "NA", "Type": "Undefined"}

        # 3 ▸ BioSample metadata
        meta = _parse_bs(_efetch("biosample", bs_acc, rettype="xml", retmode="xml"))
        host_val = meta.get("host", "")

        # 4 ▸ concatenate isolation_source-like fields
        fields = ["isolation_source", "source", "environment",
                  "sample name", "sample_name", "description"]
        iso_vals = " | ".join(meta.get(f, "") for f in fields if f in meta) or "NA"

        full_iso = f"{host_val} | {iso_vals}".strip(" |")

        return {
            "Genome": acc,
            "Isolation source": full_iso,
            "Type": classify_source(host_val, iso_vals)
        }

    except HTTPError as err:
        return {"Genome": acc,
                "Isolation source": f"ERROR: HTTP {err.code}",
                "Type": "Undefined"}
    except Exception as e:
        return {"Genome": acc,
                "Isolation source": f"ERROR: {e}",
                "Type": "Undefined"}


# ─────────────────────── main ───────────────────────
df_in = pd.read_csv(args.input_csv)
accessions = df_in[args.accession_column].dropna().astype(str).unique()

results = [classify_accession(a) for a in accessions]
pd.DataFrame(results).to_csv(args.output_csv, index=False)
print(f"[✓] Saved {len(results)} classifications → {args.output_csv}")
