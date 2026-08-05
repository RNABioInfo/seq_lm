"""Tests for the NCBI prokaryotic GTF to Oarfish converter."""

from io import StringIO
import unittest

from workflow_glue.oarfish_gtf import convert_lines, parse_attributes


INPUT_GTF = """#gtf-version 2.2
#!genome-build test
chr1\tGenbank\tgene\t10\t100\t.\t-\t.\tgene_id "coding_1"; transcript_id ""; gene "abc"; gene_biotype "protein_coding"; locus_tag "LOC1";
chr1\tGenbank\tCDS\t13\t80\t.\t-\t0\tgene_id "coding_1"; transcript_id "tx1"; protein_id "PROT1"; product "example; product";
chr1\tGenbank\tCDS\t90\t100\t.\t-\t0\tgene_id "coding_1"; transcript_id "tx1"; protein_id "PROT1";
chr1\tGenbank\tgene\t200\t250\t.\t+\t.\tgene_id "rna_1"; transcript_id ""; gene_biotype "tRNA"; locus_tag "RNA1";
chr1\tGenbank\ttranscript\t200\t250\t.\t+\t.\tgene_id "rna_1"; transcript_id "rna_tx1"; transcript_biotype "tRNA";
chr1\tGenbank\texon\t200\t250\t.\t+\t.\tgene_id "rna_1"; transcript_id "rna_tx1"; exon_number "1";
"""


class OarfishGtfConversionTest(unittest.TestCase):
    def test_conversion_synthesizes_coding_model_and_preserves_rna_model(self):
        output = StringIO()
        summary = convert_lines(INPUT_GTF.splitlines(keepends=True), output)
        converted = output.getvalue()

        self.assertEqual(summary.protein_coding_genes, 1)
        self.assertEqual(summary.synthesized_transcripts, 1)
        self.assertEqual(summary.existing_transcripts, 1)
        self.assertEqual(summary.total_transcripts, 2)
        self.assertEqual(summary.total_exons, 2)
        self.assertIn(
            'chr1\tGenbank\ttranscript\t10\t100\t.\t-\t.\tgene_id "coding_1"; '
            'transcript_id "tx1"; gene "abc"; locus_tag "LOC1"; '
            'protein_id "PROT1"; transcript_biotype "protein_coding";',
            converted,
        )
        self.assertIn(
            'chr1\tGenbank\texon\t10\t100\t.\t-\t.\tgene_id "coding_1"; '
            'transcript_id "tx1"; gene "abc"; locus_tag "LOC1"; '
            'protein_id "PROT1"; transcript_biotype "protein_coding"; '
            'exon_number "1";',
            converted,
        )
        self.assertEqual(converted.count('transcript_id "rna_tx1"'), 2)
        self.assertIn("example; product", converted)

    def test_parse_attributes_does_not_split_semicolon_inside_value(self):
        attributes = parse_attributes(
            'gene_id "g1"; note "alpha; beta"; transcript_id "tx1";'
        )
        self.assertEqual(
            attributes,
            {
                "gene_id": "g1",
                "note": "alpha; beta",
                "transcript_id": "tx1",
            },
        )

    def test_protein_coding_gene_requires_exactly_one_cds_transcript(self):
        gtf = (
            'chr1\tGenbank\tgene\t1\t20\t.\t+\t.\tgene_id "g1"; '
            'gene_biotype "protein_coding";\n'
        )
        with self.assertRaisesRegex(ValueError, "exactly one CDS transcript_id"):
            convert_lines(gtf.splitlines(keepends=True), StringIO())

    def test_existing_coding_transcript_is_not_duplicated(self):
        gtf = (
            'chr1\tGenbank\tgene\t1\t20\t.\t+\t.\tgene_id "g1"; '
            'gene_biotype "protein_coding";\n'
            'chr1\tGenbank\ttranscript\t1\t20\t.\t+\t.\tgene_id "g1"; '
            'transcript_id "tx1";\n'
            'chr1\tGenbank\tCDS\t1\t17\t.\t+\t0\tgene_id "g1"; '
            'transcript_id "tx1";\n'
        )
        with self.assertRaisesRegex(ValueError, "already has an explicit transcript"):
            convert_lines(gtf.splitlines(keepends=True), StringIO())


if __name__ == "__main__":
    unittest.main()
