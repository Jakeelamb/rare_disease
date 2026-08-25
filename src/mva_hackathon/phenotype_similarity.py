"""Local, explainable HPO-to-gene similarity using a pinned public release."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .phenotype import PhenotypeManifest, validate_phenotype_curation


@dataclass
class HpoOntology:
    parents: dict[str, set[str]]
    aliases: dict[str, str] = field(default_factory=dict)
    _cache: dict[str, frozenset[str]] = field(default_factory=dict)

    def canonical(self, term: str) -> str:
        return self.aliases.get(term, term)

    def ancestors(self, term: str, active: frozenset[str] = frozenset()) -> frozenset[str]:
        term = self.canonical(term)
        if term in self._cache:
            return self._cache[term]
        if term in active:
            raise ValueError(f"cycle in HPO ancestry at {term}")
        found = {term}
        for parent in self.parents.get(term, set()):
            found.update(self.ancestors(parent, active | {term}))
        result = frozenset(found)
        self._cache[term] = result
        return result


def load_ontology(path: Path) -> HpoOntology:
    parents: dict[str, set[str]] = {}
    aliases: dict[str, str] = {}
    current_id: str | None = None
    current_parents: set[str] = set()
    current_aliases: list[str] = []

    def finish_term() -> None:
        nonlocal current_id, current_parents, current_aliases
        if current_id:
            parents[current_id] = set(current_parents)
            aliases.update({alias: current_id for alias in current_aliases})
        current_id = None
        current_parents = set()
        current_aliases = []

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line == "[Term]":
                finish_term()
            elif line.startswith("id: HP:"):
                current_id = line.removeprefix("id: ")
            elif line.startswith("alt_id: HP:"):
                current_aliases.append(line.removeprefix("alt_id: "))
            elif line.startswith("is_a: HP:"):
                current_parents.add(line.removeprefix("is_a: ").split()[0])
            elif line.startswith("[") and line != "[Term]":
                finish_term()
        finish_term()
    return HpoOntology(parents=parents, aliases=aliases)


def load_gene_annotations(path: Path, ontology: HpoOntology) -> dict[str, set[str]]:
    annotations: dict[str, set[str]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_symbol", "hpo_id"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ValueError("gene annotation table lacks gene_symbol/hpo_id columns")
        for row in reader:
            symbol = row["gene_symbol"].strip()
            term = ontology.canonical(row["hpo_id"].strip())
            if symbol and term in ontology.parents:
                annotations.setdefault(symbol, set()).add(term)
    return annotations


def information_content(
    annotations: dict[str, set[str]], ontology: HpoOntology
) -> dict[str, float]:
    term_genes: dict[str, int] = {}
    for terms in annotations.values():
        propagated: set[str] = set()
        for term in terms:
            propagated.update(ontology.ancestors(term))
        for term in propagated:
            term_genes[term] = term_genes.get(term, 0) + 1
    total = len(annotations)
    return {term: -math.log((count + 1) / (total + 1)) for term, count in term_genes.items()}


def score_genes(
    query_terms: set[str], annotations: dict[str, set[str]], ontology: HpoOntology
) -> dict[str, float]:
    canonical_query = {ontology.canonical(term) for term in query_terms}
    missing = canonical_query - ontology.parents.keys()
    if missing:
        raise ValueError(f"HPO terms absent from pinned ontology: {sorted(missing)}")
    ic = information_content(annotations, ontology)
    query_profiles = [
        (term, ontology.ancestors(term), ic.get(term, 0.0))
        for term in sorted(canonical_query)
        if ic.get(term, 0.0) > 0
    ]
    if not query_profiles:
        raise ValueError("query terms have no informative content in gene annotations")

    term_cache: dict[str, frozenset[str]] = {}
    scores: dict[str, float] = {}
    for symbol, gene_terms in annotations.items():
        per_query: list[float] = []
        for _, query_ancestors, query_ic in query_profiles:
            best = 0.0
            for gene_term in gene_terms:
                gene_ancestors = term_cache.setdefault(gene_term, ontology.ancestors(gene_term))
                common = query_ancestors & gene_ancestors
                mica = max((ic.get(term, 0.0) for term in common), default=0.0)
                best = max(best, mica / query_ic)
            per_query.append(best)
        scores[symbol] = sum(per_query) / len(per_query)
    return scores


def build_gene_prior_document(
    phenotype: PhenotypeManifest,
    ontology_path: Path,
    annotation_path: Path,
    hpo_release: str,
) -> dict[str, Any]:
    validate_phenotype_curation(phenotype)
    ontology = load_ontology(ontology_path)
    annotations = load_gene_annotations(annotation_path, ontology)
    query = {
        observation.hpo_id
        for observation in phenotype.observations
        if observation.present and observation.subject == "proband"
    }
    scores = score_genes(query, annotations, ontology)
    return {
        "version": f"hpo-{hpo_release}",
        "assembly": "GRCh38",
        "scope": "Mechanism-blind local query-directed Resnik coverage",
        "sources": {
            "hpo": hpo_release,
            "phenotype_source_sha256": phenotype.source_sha256,
        },
        "genes": [
            {
                "symbol": symbol,
                "phenotype_gene_score": round(scores[symbol], 6),
            }
            for symbol in sorted(
                scores,
                key=lambda gene: (-scores[gene], gene),
            )
        ],
    }


def write_gene_prior_document(document: dict[str, Any], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
