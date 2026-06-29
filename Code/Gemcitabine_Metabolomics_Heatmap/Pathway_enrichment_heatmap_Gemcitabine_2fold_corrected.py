#!/usr/bin/env python3
"""
Corrected curated pathway enrichment heatmap for gemcitabine 2-fold response sets.

This version corrects the prior raw-substring annotation problem:
    if any(k in metabolite_name for k in keywords)

That previous strategy allowed short tokens such as sam, sah, amp, cmp, nad, etc.
to match inside unrelated metabolite names. For example, "sam" can be found inside
"glucosamine", and "sah" can occur as part of unrelated strings.

Correction:
1. Use explicit metabolite/component-level annotation for known metabolites.
2. Split semicolon-separated multi-component annotations and annotate each component.
3. Use token-aware regular expressions only for complete biochemical abbreviations.
4. Exclude or reassign known problematic conjugated metabolites:
   - CDP-choline/citicoline and CDP-ethanolamine -> Kennedy/choline-phospholipid, not core pyrimidine.
   - UDP-sugars -> carbohydrate/glycosylation-related class, not core pyrimidine.
   - N-acetyl-glucosamine features -> carbohydrate/amino sugar, not folate/methylation.
   - Docosahexaenoic/docosapentaenoic/stearoylethanolamide -> lipid, not folate or carnitine.
5. Export a raw-vs-corrected annotation change report for auditing.
"""

import argparse
import os
import re
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind, hypergeom
from statsmodels.stats.multitest import multipletests
from statsmodels.formula.api import ols
import statsmodels.api as sm

ALPHA = 0.05
LOG2FC_CUTOFF = 1.0

PATHWAY_ORDER = [
    "Amino acid/nitrogen",
    "Carnitine/FAO",
    "Fatty acids/lipids",
    "Folate/one-carbon/methylation",
    "Kennedy/choline-phospholipid",
    "PPP/glycolysis/carbohydrate",
    "Purine metabolism",
    "Pyrimidine metabolism",
    "Cofactors/vitamins/redox",
    "Other",
]


def normalize_text(x: str) -> str:
    s = str(x).strip().lower()
    s = s.replace("’", "'").replace("‘", "'").replace("β", "beta")
    s = re.sub(r"\s+", " ", s)
    return s


def clean_component(x: str) -> str:
    s = normalize_text(x)
    # Remove common parenthetical synonym/abbreviation chunks for exact matching,
    # but preserve full text separately by testing aliases below.
    s = re.sub(r"\s+", " ", s).strip()
    return s


def split_components(name: str):
    return [p.strip() for p in str(name).split(";") if p.strip()]


def token_match(s: str, token: str) -> bool:
    """Match a short biochemical abbreviation as a complete token only."""
    return re.search(rf"(?<![A-Za-z0-9]){re.escape(token.lower())}(?![A-Za-z0-9])", s.lower()) is not None


# Explicit component-level annotations.
# These are intentionally conservative and are designed to avoid substring false positives.
EXACT_COMPONENT_MAP = {
    # Kennedy/choline-phospholipid
    "choline": ["Kennedy/choline-phospholipid"],
    "phosphocholine": ["Kennedy/choline-phospholipid"],
    "cytidine 5'-diphosphocholine (citicoline)": ["Kennedy/choline-phospholipid"],
    "citicoline": ["Kennedy/choline-phospholipid"],
    "cdp-choline": ["Kennedy/choline-phospholipid"],
    "cdp-ethanolamine": ["Kennedy/choline-phospholipid"],
    "o-phosphorylethanolamine": ["Kennedy/choline-phospholipid"],
    "dipalmitoyl-phosphatidylcholine": ["Kennedy/choline-phospholipid", "Fatty acids/lipids"],

    # Purine metabolism
    "5'-methylthioadenosine": ["Purine metabolism"],
    "7-methylguanine": ["Purine metabolism"],
    "adenine": ["Purine metabolism"],
    "adenosine": ["Purine metabolism"],
    "adenosine 2',3'-cyclic monophosphate (2',3'-camp)": ["Purine metabolism"],
    "adenosine 3',5'-cyclic monophosphate": ["Purine metabolism"],
    "cyclic amp": ["Purine metabolism"],
    "adenosine 5'-diphosphate": ["Purine metabolism"],
    "adenosine 5'-monophosphate": ["Purine metabolism"],
    "adenosine 5'-triphosphate (atp)": ["Purine metabolism"],
    "adenosine triphosphate": ["Purine metabolism"],
    "cyclic adp-ribose": ["Purine metabolism"],
    "damp": ["Purine metabolism"],
    "dgdp": ["Purine metabolism"],
    "gtp": ["Purine metabolism"],
    "guanine": ["Purine metabolism"],
    "guanosine": ["Purine metabolism"],
    "guanosine 5'-diphosphate": ["Purine metabolism"],
    "guanosine 5'-monophosphate": ["Purine metabolism"],
    "guanosine diphosphate mannose": ["Purine metabolism", "PPP/glycolysis/carbohydrate"],
    "hypoxanthine": ["Purine metabolism"],
    "inosine": ["Purine metabolism"],
    "inosine 5'-diphosphate": ["Purine metabolism"],
    "inosine 5'-phosphate": ["Purine metabolism"],
    "xanthine": ["Purine metabolism"],
    "xanthosine": ["Purine metabolism"],
    "xanthosine 5'-monophosphate": ["Purine metabolism"],
    "1-methyladenosine": ["Purine metabolism"],
    "s-(5'-adenosyl)-l-homocysteine": ["Folate/one-carbon/methylation", "Purine metabolism"],
    "s-adenosyl methionine (sam)": ["Folate/one-carbon/methylation", "Purine metabolism"],

    # Pyrimidine metabolism
    "2'-deoxycytidine 5'-monophosphate": ["Pyrimidine metabolism"],
    "5'-cmp": ["Pyrimidine metabolism"],
    "cytidine": ["Pyrimidine metabolism"],
    "cytidine 5'-diphosphate": ["Pyrimidine metabolism"],
    "cytidine 5'-triphosphate (ctp)": ["Pyrimidine metabolism"],
    "cytidine triphosphate (ctp)": ["Pyrimidine metabolism"],
    "cytosine": ["Pyrimidine metabolism"],
    "deoxycytidine": ["Pyrimidine metabolism"],
    "orotate": ["Pyrimidine metabolism"],
    "(s)-dihydroorotate": ["Pyrimidine metabolism"],
    "dihydroorotate": ["Pyrimidine metabolism"],
    "pseudouridine": ["Pyrimidine metabolism"],
    "ribothymidine": ["Pyrimidine metabolism"],
    "thymidine": ["Pyrimidine metabolism"],
    "thymidine 5'-monophosphate": ["Pyrimidine metabolism"],
    "uracil": ["Pyrimidine metabolism"],
    "uridine": ["Pyrimidine metabolism"],
    "uridine 5'-diphosphate": ["Pyrimidine metabolism"],
    "uridine 5'-monophosphate": ["Pyrimidine metabolism"],
    "uridine triphosphate": ["Pyrimidine metabolism"],

    # UDP-sugars are treated as carbohydrate/glycosylation-associated for this class analysis.
    "uridine diphosphate glucose": ["PPP/glycolysis/carbohydrate"],
    "uridine 5'-diphosphogalactose": ["PPP/glycolysis/carbohydrate"],
    "uridine 5'-diphosphoglucuronic acid": ["PPP/glycolysis/carbohydrate"],
    "uridine 5'-diphospho-n-acetylgalactosamine": ["PPP/glycolysis/carbohydrate"],
    "uridine 5'-diphospho-n-acetylglucosamine": ["PPP/glycolysis/carbohydrate"],

    # Folate / one-carbon / methylation
    "5-methyltetrahydrofolic acid": ["Folate/one-carbon/methylation"],
    "folic acid": ["Folate/one-carbon/methylation"],
    "l-methionine": ["Folate/one-carbon/methylation", "Amino acid/nitrogen"],
    "proline betaine": ["Folate/one-carbon/methylation"],
    # Do NOT classify all N-acetyl-methionine by folate; leave as amino acid/nitrogen below.
    "n-acetyl-l-methionine": ["Amino acid/nitrogen"],

    # PPP / glycolysis / carbohydrate / amino sugar
    "6-phosphogluconic acid": ["PPP/glycolysis/carbohydrate"],
    "2-keto-3-deoxy-6-phosphogluconic acid": ["PPP/glycolysis/carbohydrate"],
    "alpha-d-galactose 1-phosphate": ["PPP/glycolysis/carbohydrate"],
    "alpha-d-glucose 1-phosphate": ["PPP/glycolysis/carbohydrate"],
    "d-fructose 6-phosphate": ["PPP/glycolysis/carbohydrate"],
    "alpha-d-galactose 1-phosphate": ["PPP/glycolysis/carbohydrate"],
    "d-(+)-trehalose": ["PPP/glycolysis/carbohydrate"],
    "d-lactose": ["PPP/glycolysis/carbohydrate"],
    "melibiose": ["PPP/glycolysis/carbohydrate"],
    "maltose": ["PPP/glycolysis/carbohydrate"],
    "d-(+)-cellobiose": ["PPP/glycolysis/carbohydrate"],
    "isomaltose": ["PPP/glycolysis/carbohydrate"],
    "d-ribose": ["PPP/glycolysis/carbohydrate"],
    "l-ribulose": ["PPP/glycolysis/carbohydrate"],
    "d-sedoheptulose": ["PPP/glycolysis/carbohydrate"],
    "fructose 1,6-bisphosphate": ["PPP/glycolysis/carbohydrate"],
    "galactitol": ["PPP/glycolysis/carbohydrate"],
    "mannitol": ["PPP/glycolysis/carbohydrate"],
    "d-sorbitol": ["PPP/glycolysis/carbohydrate"],
    "glucose 6-phosphate": ["PPP/glycolysis/carbohydrate"],
    "homovanillate": ["PPP/glycolysis/carbohydrate"],
    "3-(4-hydroxyphenyl)lactate": ["PPP/glycolysis/carbohydrate"],
    "melibiose": ["PPP/glycolysis/carbohydrate"],
    "n-acetyl-d-glucosamine": ["PPP/glycolysis/carbohydrate"],
    "n-acetyl-d-galactosamine": ["PPP/glycolysis/carbohydrate"],
    "n-acetyl-d-mannosamine": ["PPP/glycolysis/carbohydrate"],
    "n-acetylglucosamine": ["PPP/glycolysis/carbohydrate"],
    "n-acetyl glucosaminitol": ["PPP/glycolysis/carbohydrate"],
    "n-acetyl-glucosamine-6-phosphate": ["PPP/glycolysis/carbohydrate"],
    "n-acetyl-b-glucosaminylamine": ["PPP/glycolysis/carbohydrate"],
    "phosphoenolpyruvic acid": ["PPP/glycolysis/carbohydrate"],
    "ribose 1-phosphate": ["PPP/glycolysis/carbohydrate"],
    "xylulose 5-phosphate": ["PPP/glycolysis/carbohydrate"],
    "galactarate": ["PPP/glycolysis/carbohydrate"],
    "d-saccharic acid": ["PPP/glycolysis/carbohydrate"],
    "l-sorbose": ["PPP/glycolysis/carbohydrate"],
    "erytronic acid": ["PPP/glycolysis/carbohydrate"],
    "stachyose": ["PPP/glycolysis/carbohydrate"],
    "maltotetraose": ["PPP/glycolysis/carbohydrate"],
    "d-(+)-raffinose": ["PPP/glycolysis/carbohydrate"],
    "maltotriose": ["PPP/glycolysis/carbohydrate"],

    # Carnitine / FAO
    "butyryl-l-carnitine": ["Carnitine/FAO"],
    "isobutyryl-l-carnitine": ["Carnitine/FAO"],
    "decanoyl-l-carnitine": ["Carnitine/FAO"],
    "deoxycarnitine": ["Carnitine/FAO"],
    "hexanoyl-l-carnitine": ["Carnitine/FAO"],
    "isovalerylcarnitine": ["Carnitine/FAO"],
    "l-carnitine": ["Carnitine/FAO"],
    "lauroylcarnitine": ["Carnitine/FAO"],
    "lauroyl-l-carnitine": ["Carnitine/FAO"],
    "myristoyl-l-carnitine": ["Carnitine/FAO"],
    "o-acetyl-l-carnitine": ["Carnitine/FAO"],
    "octanoyl-l-carnitine": ["Carnitine/FAO"],
    "oleoyl-l-carnitine": ["Carnitine/FAO", "Fatty acids/lipids"],
    "palmitoylcarnitine": ["Carnitine/FAO"],
    "propionylcarnitine": ["Carnitine/FAO"],
    "stearoyl-l-carnitine": ["Carnitine/FAO"],

    # Fatty acids / lipids
    "arachidonic acid": ["Fatty acids/lipids"],
    "caprylic acid": ["Fatty acids/lipids"],
    "docosadienoate (22:2n6)": ["Fatty acids/lipids"],
    "docosahexaenoic acid": ["Fatty acids/lipids"],
    "docosapentaenoic acid": ["Fatty acids/lipids"],
    "eicosapentaenoic acid": ["Fatty acids/lipids"],
    "heptadecanoate": ["Fatty acids/lipids"],
    "15-methylpalmitate": ["Fatty acids/lipids"],
    "myristic acid": ["Fatty acids/lipids"],
    "palmitate": ["Fatty acids/lipids"],
    "sn-glycerol 3-phosphate": ["Fatty acids/lipids"],
    "glycerol": ["Fatty acids/lipids"],
    "glycerol 2-phosphate": ["Fatty acids/lipids"],
    "rac-glycerol 1-myristate": ["Fatty acids/lipids"],
    "stearidonic acid": ["Fatty acids/lipids"],
    "stearoylethanolamide": ["Fatty acids/lipids"],
    "adrenic acid": ["Fatty acids/lipids"],
    "10-hydroxydecanoate": ["Fatty acids/lipids"],

    # Amino acid / nitrogen
    "aspartate": ["Amino acid/nitrogen"],
    "beta-alanine": ["Amino acid/nitrogen"],
    "citrulline": ["Amino acid/nitrogen"],
    "d-tryptophan": ["Amino acid/nitrogen"],
    "guanidinoacetate": ["Amino acid/nitrogen"],
    "hypotaurine": ["Amino acid/nitrogen"],
    "l-alanine": ["Amino acid/nitrogen"],
    "d-alanine": ["Amino acid/nitrogen"],
    "l-anserine": ["Amino acid/nitrogen"],
    "l-asparagine": ["Amino acid/nitrogen"],
    "l-glutamic acid": ["Amino acid/nitrogen"],
    "o-acetyl-l-serine": ["Amino acid/nitrogen"],
    "l-glutamine": ["Amino acid/nitrogen"],
    "l-phenylalanine": ["Amino acid/nitrogen"],
    "l-serine": ["Amino acid/nitrogen"],
    "n-acetyl-l-phenylalanine": ["Amino acid/nitrogen"],
    "3-phenylpropionylglycine": ["Amino acid/nitrogen"],
    "n-alpha-acetyl-l-asparagine": ["Amino acid/nitrogen"],
    "n-formylphenylalanine": ["Amino acid/nitrogen"],
    "n-methyl-l-glutamate": ["Amino acid/nitrogen"],
    "n-methylalanine": ["Amino acid/nitrogen"],
    "3-aminoisobutanoate": ["Amino acid/nitrogen"],
    "2-amino-2-methylpropanoate": ["Amino acid/nitrogen"],
    "o-phospho-l-serine": ["Amino acid/nitrogen"],
    "o-succinyl-l-homoserine": ["Amino acid/nitrogen"],
    "taurine": ["Amino acid/nitrogen"],
    "phenylalanylalanine": ["Amino acid/nitrogen"],
    "formyl-l-methionyl peptide": ["Amino acid/nitrogen"],
    "4-aminobutanoate": ["Amino acid/nitrogen"],
    "4-guanidinobutanoate": ["Amino acid/nitrogen"],
    "n-acetylglycine": ["Amino acid/nitrogen"],

    # Cofactors / vitamins / redox
    "coenzyme a": ["Cofactors/vitamins/redox"],
    "d-pantothenic acid": ["Cofactors/vitamins/redox"],
    "dihydrobiopterin (7,8-dihydro-l-biopterin)": ["Cofactors/vitamins/redox"],
    "nad": ["Cofactors/vitamins/redox"],
    "nadp": ["Cofactors/vitamins/redox"],
    "nadph": ["Cofactors/vitamins/redox"],
    "nicotinamide": ["Cofactors/vitamins/redox"],
    "pterin": ["Cofactors/vitamins/redox"],
    "riboflavin": ["Cofactors/vitamins/redox"],
}


def corrected_pathway_annotation(metabolite_name: str) -> str:
    components = split_components(metabolite_name)
    hits = []

    for component in components:
        c = clean_component(component)

        # Exact full-component match is primary.
        if c in EXACT_COMPONENT_MAP:
            hits.extend(EXACT_COMPONENT_MAP[c])
            continue

        # Some features contain a primary name with a parenthetical alias.
        c_no_paren = re.sub(r"\s*\([^)]*\)\s*", " ", c)
        c_no_paren = re.sub(r"\s+", " ", c_no_paren).strip()
        if c_no_paren in EXACT_COMPONENT_MAP:
            hits.extend(EXACT_COMPONENT_MAP[c_no_paren])
            continue

        # Conservative fallback for complete abbreviation tokens only.
        # These are not raw substrings.
        for tok in ["amp", "adp", "atp", "gmp", "gdp", "gtp", "imp", "idp", "xmp", "damp", "dgmp", "dgdp"]:
            if token_match(c, tok):
                hits.append("Purine metabolism")
                break

        for tok in ["ump", "udp", "utp", "cmp", "cdp", "ctp", "tmp", "dtmp", "dcmp", "dfdcmp"]:
            if token_match(c, tok):
                # Special cases: Kennedy CDP-conjugates and UDP-sugars should not be core pyrimidine calls.
                if any(term in c for term in ["choline", "ethanolamine"]):
                    hits.append("Kennedy/choline-phospholipid")
                elif any(term in c for term in ["glucose", "galactose", "glucuronic", "acetylglucosamine", "acetylgalactosamine", "sugar"]):
                    hits.append("PPP/glycolysis/carbohydrate")
                else:
                    hits.append("Pyrimidine metabolism")
                break

        # Conservative long-token fallback. No short tokens like sam/sah/nad/amp by raw substring.
        if re.search(r"\b(orotate|dihydroorotate|uridine|cytidine|cytosine|uracil|thymidine|pseudouridine|ribothymidine)\b", c):
            hits.append("Pyrimidine metabolism")
        if re.search(r"\b(adenosine|adenine|guanine|guanosine|inosine|hypoxanthine|xanthine|xanthosine)\b", c):
            hits.append("Purine metabolism")
        if re.search(r"\b(folic acid|methyltetrahydrofolic|s-adenosyl|homocysteine|betaine)\b", c):
            hits.append("Folate/one-carbon/methylation")
        if re.search(r"\b(phosphocholine|citicoline|phosphorylethanolamine|phosphoethanolamine|phosphatidylcholine|choline)\b", c):
            hits.append("Kennedy/choline-phospholipid")
        if re.search(r"\b(glucose|fructose|phosphogluconic|sedoheptulose|galactitol|mannitol|sorbitol|trehalose|maltose|melibiose|ribose|xylulose|raffinose|stachyose|sorbose|glucosamine|galactosamine)\b", c):
            hits.append("PPP/glycolysis/carbohydrate")
        if re.search(r"\b(carnitine|carnitines|lauroylcarnitine|isovalerylcarnitine|palmitoylcarnitine|propionylcarnitine|deoxycarnitine)\b", c):
            hits.append("Carnitine/FAO")
        if re.search(r"\b(arachidonic|palmitate|palmitic|stearidonic|docosahexaenoic|docosapentaenoic|eicosapentaenoic|myristic|glycerol|docosadienoate|caprylic|adrenic|hydroxydecanoate|stearoylethanolamide)\b", c):
            hits.append("Fatty acids/lipids")
        if re.search(r"\b(glutamine|glutamate|aspartate|asparagine|tryptophan|phenylalanine|serine|alanine|citrulline|guanidino|taurine|hypotaurine|aminobutanoate|beta-alanine|methionyl|methionine|acetylglycine|phenylalanylalanine)\b", c):
            # Avoid moving SAM/SAH away from folate; keep both if exact/fallback already added.
            hits.append("Amino acid/nitrogen")
        if re.search(r"\b(pantothenic|riboflavin|nicotinamide|coenzyme|pterin|biopterin)\b", c) or token_match(c, "nad") or token_match(c, "nadp") or token_match(c, "nadph"):
            hits.append("Cofactors/vitamins/redox")

    # Remove duplicates while keeping pathway order.
    unique_hits = [p for p in PATHWAY_ORDER if p in set(hits) and p != "Other"]
    return "; ".join(unique_hits) if unique_hits else "Other"


def raw_substring_annotation(metabolite_name: str) -> str:
    """Replicate the original permissive substring logic for auditing only."""
    pathways = {
        "Kennedy/choline-phospholipid": [
            "phosphocholine", "citicoline", "cdp-choline", "diphosphocholine",
            "phosphorylethanolamine", "cdp-ethanolamine", "choline"
        ],
        "Purine metabolism": [
            "adenosine", "adenine", "amp", "adp", "atp", "guanine", "guanosine",
            "gmp", "gdp", "gtp", "inosine", "imp", "hypoxanthine", "xanthosine",
            "xanthine", "damp", "dgmp"
        ],
        "Pyrimidine metabolism": [
            "uridine", "uracil", "ump", "udp", "utp", "cytidine", "cytosine",
            "cmp", "cdp", "ctp", "thymidine", "tmp", "dtmp", "orotate",
            "dihydroorotate", "dfdcmp"
        ],
        "Folate/one-carbon/methylation": [
            "folic", "methyltetrahydrofolic", "methionine", "s-adenosyl",
            "homocysteine", "sah", "sam", "betaine"
        ],
        "PPP/glycolysis/carbohydrate": [
            "glucose", "fructose", "phosphogluconic", "sedoheptulose", "ribose",
            "lactate", "pyruvic", "phosphoenolpyruvic", "galactitol", "mannitol",
            "sorbitol", "trehalose", "maltose", "melibiose"
        ],
        "Carnitine/FAO": [
            "carnitine", "acylcarnitine", "isovaleryl", "lauroyl", "decanoyl",
            "hexanoyl", "myristoyl", "stearoyl", "propionylcarnitine"
        ],
        "Fatty acids/lipids": [
            "palmitate", "palmitic", "stearidonic", "arachidonic", "myristic",
            "oleoyl", "glycerol", "dihome", "eicosenoic", "docosadienoate",
            "caprylic", "linoleic"
        ],
        "Amino acid/nitrogen": [
            "glutamine", "glutamate", "aspartate", "asparagine", "tryptophan",
            "phenylalanine", "serine", "alanine", "citrulline", "guanidino",
            "taurine", "hypotaurine", "aminobutanoate", "beta-alanine"
        ],
        "Cofactors/vitamins/redox": [
            "nadp", "nad", "pantothenic", "riboflavin", "nicotinamide",
            "coenzyme", "pterin"
        ],
    }
    n = str(metabolite_name).lower()
    hits = []
    for pathway, keywords in pathways.items():
        if any(k.lower() in n for k in keywords):
            hits.append(pathway)
    return "; ".join(hits) if hits else "Other"


def load_and_preprocess(input_file):
    df_raw = pd.read_excel(input_file, sheet_name=0)
    id_candidates = ["row identity (all IDs)", "Metabolite ID", "Metabolite", "metabolite"]
    id_col = next((c for c in id_candidates if c in df_raw.columns), df_raw.columns[0])

    replicate_pattern = re.compile(r"^(2N|4N)_[CG][1-4]$")
    rep_cols = [c for c in df_raw.columns if replicate_pattern.match(str(c))]
    expected_cols = [f"{p}_{t}{i}" for p in ["2N", "4N"] for t in ["C", "G"] for i in range(1, 5)]
    ordered_cols = [c for c in expected_cols if c in rep_cols]

    if len(ordered_cols) != 16:
        raise ValueError(f"Expected 16 replicate columns, found {len(ordered_cols)}: {ordered_cols}")

    df = df_raw[[id_col] + ordered_cols].copy()
    df[id_col] = df[id_col].astype(str).str.strip()
    counts = df[id_col].value_counts()
    df["Feature_ID"] = [
        f"{name} | row_{i}" if counts[name] > 1 else name
        for i, name in zip(df_raw.index, df[id_col])
    ]

    X_raw = df[ordered_cols].apply(pd.to_numeric, errors="coerce")
    keep = (X_raw.fillna(0) > 0).sum(axis=1) > 0
    df = df.loc[keep].reset_index(drop=True)
    X_raw = X_raw.loc[keep].reset_index(drop=True)

    X_imp = X_raw.copy()
    row_min_positive = X_imp.where(X_imp > 0).min(axis=1)
    global_min_positive = X_imp.where(X_imp > 0).min().min()
    row_min_positive = row_min_positive.fillna(global_min_positive)

    for i in range(X_imp.shape[0]):
        fill_value = row_min_positive.iloc[i] / 2.0
        X_imp.iloc[i, :] = X_imp.iloc[i, :].replace(0, np.nan).fillna(fill_value)

    X_log2 = np.log2(X_imp)
    X_log2.index = df["Feature_ID"]

    meta = []
    for col in ordered_cols:
        ploidy, rest = col.split("_")
        treatment = "Gemcitabine" if rest.startswith("G") else "Control"
        group = f"{ploidy}_{'G' if treatment == 'Gemcitabine' else 'C'}"
        meta.append({
            "sample": col,
            "ploidy": ploidy,
            "treatment": treatment,
            "group": group,
            "replicate": rest[1:],
        })
    meta = pd.DataFrame(meta)

    group_cols = {
        "2N_C": [c for c in ordered_cols if c.startswith("2N_C")],
        "2N_G": [c for c in ordered_cols if c.startswith("2N_G")],
        "4N_C": [c for c in ordered_cols if c.startswith("4N_C")],
        "4N_G": [c for c in ordered_cols if c.startswith("4N_G")],
    }
    return df, X_log2, meta, group_cols, ordered_cols, id_col


def compute_response_sets(df, X_log2, meta, group_cols, id_col):
    log2fc_2n = X_log2[group_cols["2N_G"]].mean(axis=1) - X_log2[group_cols["2N_C"]].mean(axis=1)
    log2fc_4n = X_log2[group_cols["4N_G"]].mean(axis=1) - X_log2[group_cols["4N_C"]].mean(axis=1)

    p_2n = ttest_ind(
        X_log2[group_cols["2N_G"]].values,
        X_log2[group_cols["2N_C"]].values,
        axis=1,
        equal_var=False,
        nan_policy="omit",
    ).pvalue
    p_4n = ttest_ind(
        X_log2[group_cols["4N_G"]].values,
        X_log2[group_cols["4N_C"]].values,
        axis=1,
        equal_var=False,
        nan_policy="omit",
    ).pvalue

    p_interaction = []
    for idx in range(X_log2.shape[0]):
        dat = meta.copy()
        dat["value"] = X_log2.iloc[idx].values
        try:
            model = ols(
                "value ~ C(ploidy) + C(treatment) + C(ploidy):C(treatment)",
                data=dat,
            ).fit()
            aov = sm.stats.anova_lm(model, typ=2)
            p_interaction.append(aov.loc["C(ploidy):C(treatment)", "PR(>F)"])
        except Exception:
            p_interaction.append(np.nan)

    p_interaction = np.asarray(p_interaction, dtype=float)
    diff_response = log2fc_2n.values - log2fc_4n.values

    results = pd.DataFrame({
        "Feature_ID": df["Feature_ID"],
        "Metabolite": df[id_col],
        "log2FC_2N_G_vs_C": log2fc_2n.values,
        "p_2N_G_vs_C": p_2n,
        "log2FC_4N_G_vs_C": log2fc_4n.values,
        "p_4N_G_vs_C": p_4n,
        "Differential_response_2N_minus_4N": diff_response,
        "p_interaction": p_interaction,
    })

    sets = {
        "2N_up_after_gem_2fold": (log2fc_2n.values >= LOG2FC_CUTOFF) & (p_2n < ALPHA),
        "2N_down_after_gem_2fold": (log2fc_2n.values <= -LOG2FC_CUTOFF) & (p_2n < ALPHA),
        "4N_up_after_gem_2fold": (log2fc_4n.values >= LOG2FC_CUTOFF) & (p_4n < ALPHA),
        "4N_down_after_gem_2fold": (log2fc_4n.values <= -LOG2FC_CUTOFF) & (p_4n < ALPHA),
        "Differential_response_2fold": (np.abs(diff_response) >= LOG2FC_CUTOFF) & (p_interaction < ALPHA),
        "Interaction_p05": (p_interaction < ALPHA),
    }
    return results, sets


def expand_classes(df, class_col="Corrected_pathway_class"):
    rows = []
    for _, r in df.iterrows():
        classes = [c.strip() for c in str(r[class_col]).split(";") if c.strip()]
        if not classes:
            classes = ["Other"]
        for c in classes:
            rr = r.to_dict()
            rr["Expanded_pathway_class"] = c
            rows.append(rr)
    return pd.DataFrame(rows)


def hypergeometric_enrichment(results, sets):
    results = results.copy()
    results["Raw_substring_pathway_class"] = results["Metabolite"].apply(raw_substring_annotation)
    results["Corrected_pathway_class"] = results["Metabolite"].apply(corrected_pathway_annotation)
    results["annotation_changed"] = results["Raw_substring_pathway_class"] != results["Corrected_pathway_class"]

    # Expanded rows allow multi-class features to contribute one count to each annotated class,
    # while each feature is counted at most once per pathway.
    expanded_all = expand_classes(results, "Corrected_pathway_class")
    M = len(results)

    rows = []
    for set_name, mask in sets.items():
        selected = results.loc[mask].copy()
        selected_expanded = expand_classes(selected, "Corrected_pathway_class")
        N = len(selected)

        for pathway in PATHWAY_ORDER:
            if pathway == "Other":
                bg_features = set(expanded_all.loc[expanded_all["Expanded_pathway_class"] == "Other", "Feature_ID"])
                sel_features = set(selected_expanded.loc[selected_expanded["Expanded_pathway_class"] == "Other", "Feature_ID"])
            else:
                bg_features = set(expanded_all.loc[expanded_all["Expanded_pathway_class"] == pathway, "Feature_ID"])
                sel_features = set(selected_expanded.loc[selected_expanded["Expanded_pathway_class"] == pathway, "Feature_ID"])

            K = len(bg_features)
            x = len(sel_features)
            p_value = hypergeom.sf(x - 1, M, K, N) if x > 0 and K > 0 and N > 0 else 1.0

            rows.append({
                "Set": set_name,
                "Pathway_class": pathway,
                "Selected_count": int(x),
                "Background_count": int(K),
                "Set_size": int(N),
                "Background_size": int(M),
                "p_value": p_value,
            })

    enrichment = pd.DataFrame(rows)
    enrichment["FDR_BH"] = np.nan

    for set_name in enrichment["Set"].unique():
        idx = enrichment["Set"] == set_name
        enrichment.loc[idx, "FDR_BH"] = multipletests(
            enrichment.loc[idx, "p_value"].astype(float), method="fdr_bh"
        )[1]

    return enrichment, results


def make_heatmap(enrichment, output_png, output_pdf=None, output_matrix_csv=None):
    enr = enrichment[
        (enrichment["Pathway_class"] != "Other")
        & (enrichment["Selected_count"] > 0)
    ].copy()

    enr["score"] = -np.log10(np.clip(enr["FDR_BH"].astype(float), 1e-300, None))
    enr_show = enr[(enr["p_value"] < 0.25) | (enr["FDR_BH"] < 0.25)].copy()
    if len(enr_show) == 0:
        enr_show = enr.sort_values("p_value").head(25)

    pivot = enr_show.pivot_table(
        index="Pathway_class", columns="Set", values="score", aggfunc="max"
    ).fillna(0)

    set_order = [
        "2N_up_after_gem_2fold",
        "2N_down_after_gem_2fold",
        "4N_up_after_gem_2fold",
        "4N_down_after_gem_2fold",
        "Differential_response_2fold",
        "Interaction_p05",
    ]
    pivot = pivot[[s for s in set_order if s in pivot.columns]]

    # Biology-oriented row ordering.
    row_order = [
        "Amino acid/nitrogen",
        "Carnitine/FAO",
        "Fatty acids/lipids",
        "Folate/one-carbon/methylation",
        "Kennedy/choline-phospholipid",
        "PPP/glycolysis/carbohydrate",
        "Purine metabolism",
        "Pyrimidine metabolism",
        "Cofactors/vitamins/redox",
    ]
    pivot = pivot.loc[[r for r in row_order if r in pivot.index]]

    if output_matrix_csv:
        pivot.to_csv(output_matrix_csv)

    fig, ax = plt.subplots(figsize=(12, max(5, 0.50 * len(pivot))))
    im = ax.imshow(pivot.values, aspect="auto", cmap="viridis", interpolation="nearest")
    ax.set_title("Corrected curated pathway-class enrichment across 2-fold response sets", fontsize=14)
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=90)
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels(pivot.index)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("-log10(FDR)")
    ax.set_xlabel("Response set")
    ax.set_ylabel("Corrected pathway / metabolite class")
    plt.tight_layout()
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    if output_pdf:
        fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)
    return pivot


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input metabolomics .xlsm/.xlsx file")
    parser.add_argument("--outdir", help="Output directory")
    parser.add_argument("--output-dir", help="Canonical output directory alias for --outdir.")
    args = parser.parse_args()
    if args.output_dir:
        args.outdir = args.output_dir
    if not args.outdir:
        parser.error("--outdir or --output-dir is required")

    os.makedirs(args.outdir, exist_ok=True)
    figdir = os.path.join(args.outdir, "figures")
    tabledir = os.path.join(args.outdir, "tables")
    os.makedirs(figdir, exist_ok=True)
    os.makedirs(tabledir, exist_ok=True)

    df, X_log2, meta, group_cols, ordered_cols, id_col = load_and_preprocess(args.input)
    results, sets = compute_response_sets(df, X_log2, meta, group_cols, id_col)
    enrichment, annotated_results = hypergeometric_enrichment(results, sets)

    annotated_results.to_csv(os.path.join(tabledir, "annotated_metabolites_with_response_stats_corrected.csv"), index=False)
    enrichment.to_csv(os.path.join(tabledir, "corrected_curated_pathway_enrichment_2fold.csv"), index=False)

    changed = annotated_results[annotated_results["annotation_changed"]].copy()
    changed.to_csv(os.path.join(tabledir, "annotation_change_report_raw_substring_vs_corrected.csv"), index=False)

    # Summary tables.
    annotation_counts = (
        expand_classes(annotated_results, "Corrected_pathway_class")
        .groupby("Expanded_pathway_class")["Feature_ID"].nunique()
        .reset_index(name="n_features")
        .sort_values("n_features", ascending=False)
    )
    annotation_counts.to_csv(os.path.join(tabledir, "corrected_annotation_background_counts.csv"), index=False)

    # Changed pathways summary.
    change_summary = (
        changed.groupby(["Raw_substring_pathway_class", "Corrected_pathway_class"])
        .size()
        .reset_index(name="n_features")
        .sort_values("n_features", ascending=False)
    )
    change_summary.to_csv(os.path.join(tabledir, "annotation_change_summary.csv"), index=False)

    matrix = make_heatmap(
        enrichment,
        output_png=os.path.join(figdir, "corrected_curated_pathway_enrichment_heatmap_2fold.png"),
        output_pdf=os.path.join(figdir, "corrected_curated_pathway_enrichment_heatmap_2fold.pdf"),
        output_matrix_csv=os.path.join(tabledir, "corrected_curated_pathway_enrichment_heatmap_matrix.csv"),
    )

    readme = f"""
Corrected pathway enrichment heatmap
====================================

Main correction:
The prior script used permissive raw substring matching:
    if any(k.lower() in metabolite_name for k in keywords)

This could falsely assign pathways when short tokens such as sam, sah, amp, cmp, udp, or nad
appeared inside unrelated metabolite names. This script replaces that with explicit
component-level annotations and token-aware matching for complete biochemical abbreviations.

Input features after preprocessing: {len(annotated_results)}
Features with changed annotation compared with old raw-substring approach: {len(changed)}

Main outputs:
- figures/corrected_curated_pathway_enrichment_heatmap_2fold.png
- figures/corrected_curated_pathway_enrichment_heatmap_2fold.pdf
- tables/corrected_curated_pathway_enrichment_2fold.csv
- tables/corrected_curated_pathway_enrichment_heatmap_matrix.csv
- tables/annotated_metabolites_with_response_stats_corrected.csv
- tables/annotation_change_report_raw_substring_vs_corrected.csv
- tables/annotation_change_summary.csv

Interpretation:
This is a corrected curated, name-based pathway-class enrichment. It is still not a substitute
for complete KEGG/HMDB ID-mapped pathway enrichment, but it removes the most important false
positive mechanism from raw substring matching.
"""
    with open(os.path.join(args.outdir, "README_corrected_pathway_enrichment.txt"), "w") as f:
        f.write(readme)

    print(f"Loaded {len(annotated_results)} metabolite features.")
    print(f"Changed annotations after correction: {len(changed)}")
    print(f"Heatmap matrix shape: {matrix.shape[0]} pathway classes x {matrix.shape[1]} response sets.")
    print("Saved corrected heatmap and tables to:", args.outdir)


if __name__ == "__main__":
    main()
