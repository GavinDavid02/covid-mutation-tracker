"""
tests/test_pipeline.py
──────────────────────
Unit tests for the COVID-19 Mutation Tracker pipeline.

Run:
    python -m pytest tests/ -v
    python -m pytest tests/ -v --tb=short
"""

import json
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))


class TestConfig(unittest.TestCase):
    """Test configuration constants and paths."""

    def test_imports(self):
        from config import (
            SPIKE_START, SPIKE_END, SPIKE_DOMAINS,
            KEY_MUTATIONS, VARIANT_DEFINITIONS,
        )
        self.assertIsInstance(SPIKE_START, int)
        self.assertIsInstance(SPIKE_END,   int)
        self.assertGreater(SPIKE_END, SPIKE_START)

    def test_spike_coordinates(self):
        from config import SPIKE_START, SPIKE_END
        self.assertEqual(SPIKE_START, 21562)
        self.assertEqual(SPIKE_END,   25384)

    def test_spike_domains_coverage(self):
        from config import SPIKE_DOMAINS, SPIKE_AA_LEN
        # Every domain must be within protein bounds
        for domain, (start, end) in SPIKE_DOMAINS.items():
            self.assertGreater(start, 0, f"{domain} start must be > 0")
            self.assertLessEqual(end, SPIKE_AA_LEN,
                                 f"{domain} end {end} > spike length {SPIKE_AA_LEN}")
            self.assertLess(start, end,
                            f"{domain} start {start} >= end {end}")

    def test_key_mutations_format(self):
        from config import KEY_MUTATIONS
        for entry in KEY_MUTATIONS:
            self.assertEqual(len(entry), 4,
                             f"KEY_MUTATIONS entry should have 4 fields: {entry}")
            name, pos, cat, _ = entry
            self.assertIsInstance(name, str)
            self.assertIsInstance(pos, int)
            self.assertGreater(pos, 0)

    def test_variant_definitions(self):
        from config import VARIANT_DEFINITIONS
        required_keys = {"pango", "who", "defining_spike", "origin_country",
                         "first_detected", "colour"}
        for var_name, v in VARIANT_DEFINITIONS.items():
            missing = required_keys - set(v.keys())
            self.assertFalse(missing,
                             f"Variant {var_name} missing keys: {missing}")
            self.assertTrue(v["colour"].startswith("#"),
                            f"Colour {v['colour']} should be a hex string")


class TestMutationParsing(unittest.TestCase):
    """Test amino-acid substitution parsing from Nextstrain metadata."""

    def setUp(self):
        from pipeline.alignment import parse_aa_substitutions, parse_aa_deletions
        self.parse_subs = parse_aa_substitutions
        self.parse_dels = parse_aa_deletions

    def test_parse_single_mutation(self):
        result = self.parse_subs("S:D614G", gene="S")
        self.assertEqual(result, ["D614G"])

    def test_parse_multiple_mutations(self):
        result = self.parse_subs("S:N501Y,S:P681H,S:D614G,ORF1a:T4715I", gene="S")
        self.assertIn("N501Y", result)
        self.assertIn("P681H", result)
        self.assertIn("D614G", result)
        self.assertNotIn("T4715I", result)   # ORF1a, not spike

    def test_parse_empty_string(self):
        self.assertEqual(self.parse_subs(""),   [])
        self.assertEqual(self.parse_subs(None), [])

    def test_parse_deletions(self):
        result = self.parse_dels("S:HV69-70,ORF1a:3675-3677", gene="S")
        self.assertIn("HV69-70del", result)
        self.assertEqual(len(result), 1)   # ORF1a deletion excluded


class TestPositionHelpers(unittest.TestCase):
    """Test mutation string position extraction and domain lookup."""

    def test_extract_position(self):
        from pipeline.alignment import _extract_position
        self.assertEqual(_extract_position("D614G"),  614)
        self.assertEqual(_extract_position("N501Y"),  501)
        self.assertEqual(_extract_position("K417N"),  417)
        self.assertIsNone(_extract_position("delHV"))

    def test_position_to_domain(self):
        from pipeline.alignment import _position_to_domain
        # RBD: 306-527
        self.assertEqual(_position_to_domain(417), "RBD")
        self.assertEqual(_position_to_domain(501), "RBD")
        # NTD: 14-305
        self.assertEqual(_position_to_domain(144), "NTD")
        # Signal peptide: 1-13
        self.assertEqual(_position_to_domain(5),   "Signal peptide")
        # Beyond protein
        self.assertEqual(_position_to_domain(2000), "Other")

    def test_lineage_to_who(self):
        from pipeline.alignment import _lineage_to_who
        self.assertEqual(_lineage_to_who("B.1.1.7"),   "Alpha")
        self.assertEqual(_lineage_to_who("B.1.617.2"), "Delta")
        self.assertEqual(_lineage_to_who("BA.1"),      "Omicron")
        self.assertEqual(_lineage_to_who("BA.2"),      "Omicron")
        self.assertEqual(_lineage_to_who("XBB.1.5"),   "Omicron")
        self.assertEqual(_lineage_to_who("JN.1"),      "Omicron")
        self.assertEqual(_lineage_to_who("Unknown"),   "Ancestral")


class TestFrequencyAnalysis(unittest.TestCase):
    """Test frequency computation functions."""

    def _make_df(self, n=200):
        """Build a minimal synthetic mutation DataFrame for testing."""
        import random
        random.seed(0)
        np.random.seed(0)

        dates = pd.date_range("2021-01-01", "2022-12-31", periods=n)
        labels = (["Alpha"] * (n // 4) + ["Delta"] * (n // 4) +
                  ["Omicron"] * (n // 2))
        random.shuffle(labels)

        records = []
        for i in range(n):
            who = labels[i]
            if who == "Alpha":
                muts = ["D614G", "N501Y", "P681H"]
            elif who == "Delta":
                muts = ["D614G", "L452R", "T478K"]
            else:
                muts = ["D614G", "K417N", "E484A", "N501Y"]

            records.append({
                "strain":          f"seq_{i}",
                "date":            dates[i],
                "region":          "Europe",
                "country":         "United Kingdom",
                "lineage":         who.lower(),
                "who_label":       who,
                "spike_mutations": json.dumps(muts),
                "spike_mut_count": len(muts),
                "rbd_mutations":   json.dumps([m for m in muts if "N501" in m or "K417" in m]),
                "rbd_mut_count":   1,
                **{f"has_{m}": (m in muts) for m in ["D614G", "N501Y", "L452R", "K417N", "E484A", "T478K", "P681H"]},
            })

        df = pd.DataFrame(records)
        df["date"] = pd.to_datetime(df["date"])
        df["year_month"] = df["date"].dt.to_period("M")
        return df

    def test_variant_frequency_shape(self):
        from analysis.frequency_analysis import compute_variant_frequency
        df = self._make_df()
        result = compute_variant_frequency(df)
        self.assertIn("year_month",  result.columns)
        self.assertIn("who_label",   result.columns)
        self.assertIn("frequency",   result.columns)
        self.assertGreater(len(result), 0)

    def test_frequency_between_0_and_100(self):
        from analysis.frequency_analysis import compute_variant_frequency
        df = self._make_df()
        result = compute_variant_frequency(df)
        self.assertTrue((result["frequency"] >= 0).all())
        self.assertTrue((result["frequency"] <= 100).all())

    def test_mutation_frequency_time(self):
        from analysis.frequency_analysis import compute_mutation_frequency_over_time
        df = self._make_df()
        result = compute_mutation_frequency_over_time(df, mutations=["D614G", "N501Y"])
        self.assertIn("mutation", result.columns)
        self.assertIn("frequency", result.columns)
        # D614G should be near 100% (present in all synthetic sequences)
        d614g = result[result["mutation"] == "D614G"]["frequency"]
        self.assertGreater(d614g.mean(), 80)

    def test_convergent_evolution(self):
        from analysis.frequency_analysis import detect_convergent_evolution
        df = self._make_df()
        result = detect_convergent_evolution(df)
        self.assertIn("mutation",         result.columns)
        self.assertIn("appearance_count", result.columns)
        # D614G should appear in all variants = convergent
        if not result.empty:
            top = result.iloc[0]["mutation"]
            self.assertIn("D614G", result["mutation"].values)

    def test_hotspot_detection(self):
        from analysis.frequency_analysis import compute_mutation_hotspots
        df = self._make_df()
        result = compute_mutation_hotspots(df)
        self.assertIn("position",       result.columns)
        self.assertIn("mutation_count", result.columns)
        self.assertIn("is_hotspot",     result.columns)
        # Position 614 (D614G) should be highest
        top_pos = result.nlargest(1, "mutation_count")["position"].iloc[0]
        self.assertEqual(top_pos, 614)

    def test_cooccurrence(self):
        from analysis.frequency_analysis import compute_mutation_cooccurrence
        df = self._make_df(300)
        result = compute_mutation_cooccurrence(df, top_n=10)
        if not result.empty:
            self.assertIn("mut_a", result.columns)
            self.assertIn("mut_b", result.columns)
            self.assertIn("odds_ratio", result.columns)
            # All odds ratios should be positive
            self.assertTrue((result["odds_ratio"] > 0).all())


class TestSampleDataGenerator(unittest.TestCase):
    """Test synthetic data generation."""

    def test_generate_sample_metadata(self):
        """Generated data should have the right shape and columns."""
        from pipeline.fetch_data import _generate_sample_metadata
        from config import OUTPUT_DIR
        import os, tempfile

        # Generate small sample
        orig_dir = OUTPUT_DIR
        import config as cfg
        tmp = Path(tempfile.mkdtemp())
        cfg.OUTPUT_DIR = tmp  # redirect output
        try:
            path = _generate_sample_metadata(n=200)
            df = pd.read_parquet(path)
            self.assertGreater(len(df), 100)
            for col in ["strain", "date", "country", "pangolin_lineage",
                        "aaSubstitutions"]:
                self.assertIn(col, df.columns, f"Missing column: {col}")
            # Dates should be parseable
            dates = pd.to_datetime(df["date"])
            self.assertFalse(dates.isna().all())
        finally:
            cfg.OUTPUT_DIR = orig_dir

    def test_lineage_variety(self):
        """Generated data should include multiple lineages."""
        from pipeline.fetch_data import _generate_sample_metadata
        from config import OUTPUT_DIR
        import config as cfg
        import tempfile

        tmp = Path(tempfile.mkdtemp())
        cfg.OUTPUT_DIR = tmp
        try:
            path = _generate_sample_metadata(n=1000)
            df = pd.read_parquet(path)
            lineages = df["pangolin_lineage"].unique()
            self.assertGreater(len(lineages), 5, "Expected >5 lineages in 1000 samples")
        finally:
            cfg.OUTPUT_DIR = OUTPUT_DIR


class TestVisualization(unittest.TestCase):
    """Test that visualization functions return valid Plotly figures."""

    def _sample_variant_freq(self):
        months = pd.period_range("2021-01", "2022-12", freq="M")
        rows = []
        for m in months:
            for who, freq in [("Alpha", 40), ("Delta", 35), ("Omicron", 25)]:
                rows.append({"year_month": m, "who_label": who,
                             "count": 10, "total": 30, "frequency": freq})
        return pd.DataFrame(rows)

    def test_timeline_variant(self):
        from visualization.plots import timeline_plot
        vf = self._sample_variant_freq()
        fig = timeline_plot(vf, mode="variant")
        self.assertIsNotNone(fig)
        self.assertGreater(len(fig.data), 0)

    def test_phylogenetic_tree(self):
        from visualization.plots import phylogenetic_tree_plotly
        fig = phylogenetic_tree_plotly()
        self.assertIsNotNone(fig)
        self.assertGreater(len(fig.data), 0)

    def test_spike_protein_map_empty(self):
        """Spike map should not crash on empty hotspot data."""
        from visualization.plots import spike_protein_map
        fig = spike_protein_map(pd.DataFrame())
        self.assertIsNotNone(fig)

    def test_spike_protein_map_data(self):
        """Spike map with real-ish hotspot data."""
        from visualization.plots import spike_protein_map
        hs = pd.DataFrame({
            "position":       [417, 484, 501, 614, 681],
            "mutation_count": [50, 70, 80, 200, 60],
            "top_mutations":  ["K417N", "E484K", "N501Y", "D614G", "P681H"],
            "domain":         ["RBD", "RBD", "RBD", "Central domain", "Furin cleavage"],
            "is_hotspot":     [False, True, True, True, False],
        })
        fig = spike_protein_map(hs)
        self.assertIsNotNone(fig)
        self.assertGreater(len(fig.data), 0)


class TestEndToEnd(unittest.TestCase):
    """
    Lightweight end-to-end smoke test: generate data → align → analyse.
    Uses only 100 sequences to keep it fast.
    """

    def test_full_pipeline_smoke(self):
        import tempfile, config as cfg
        from pathlib import Path

        tmp = Path(tempfile.mkdtemp())
        cfg.OUTPUT_DIR = tmp
        cfg.RAW_DIR    = tmp / "raw"
        cfg.RAW_DIR.mkdir()

        try:
            from pipeline.fetch_data import _generate_sample_metadata
            meta_path = _generate_sample_metadata(n=100)

            from pipeline.alignment import build_mutation_table
            df = build_mutation_table(meta_path)
            self.assertGreater(len(df), 50)

            from analysis.frequency_analysis import run_all_analyses
            results = run_all_analyses(df)
            self.assertIn("variant_frequency", results)
            self.assertIn("convergent", results)
            self.assertIn("hotspots", results)

        finally:
            from config import OUTPUT_DIR as real_out, RAW_DIR as real_raw
            cfg.OUTPUT_DIR = real_out
            cfg.RAW_DIR    = real_raw


if __name__ == "__main__":
    unittest.main(verbosity=2)
