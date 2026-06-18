#!/usr/bin/env python3
"""Check assembly ends for Streptomycete telomere sequences using BLAST.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


@dataclass(frozen=True)
class FastaRecord:
	name: str
	sequence: str


@dataclass(frozen=True)
class BlastHit:
	query_id: str
	subject_id: str
	pident: float
	length: int
	qstart: int
	qend: int
	sstart: int
	send: int
	evalue: float
	bitscore: float
	qlen: int
	slen: int

	@property
	def qcov_pct(self) -> float:
		return 100.0 * self.length / self.qlen if self.qlen else 0.0

	@property
	def strand(self) -> str:
		return "+" if self.send >= self.sstart else "-"

	@property
	def subject_start(self) -> int:
		return min(self.sstart, self.send)

	@property
	def subject_end(self) -> int:
		return max(self.sstart, self.send)

	def subject_covers_start(self, length: int) -> bool:
		return self.subject_start <= 1 and self.subject_end >= length


def read_fasta(path: Path) -> Iterator[FastaRecord]:
	name = None
	chunks: list[str] = []
	with path.open() as handle:
		for raw_line in handle:
			line = raw_line.strip()
			if not line:
				continue
			if line.startswith(">"):
				if name is not None:
					yield FastaRecord(name=name, sequence="".join(chunks))
				name = line[1:].split(None, 1)[0]
				chunks = []
			else:
				chunks.append(line)
	if name is not None:
		yield FastaRecord(name=name, sequence="".join(chunks))


def write_fasta(records: Iterable[FastaRecord], path: Path) -> None:
	with path.open("w") as handle:
		for record in records:
			handle.write(f">{record.name}\n")
			sequence = record.sequence
			for offset in range(0, len(sequence), 60):
				handle.write(sequence[offset : offset + 60] + "\n")


def build_blast_db(telomere_fasta: Path, db_prefix: Path) -> None:
	command = [
		"makeblastdb",
		"-in",
		str(telomere_fasta),
		"-dbtype",
		"nucl",
		"-parse_seqids",
		"-out",
		str(db_prefix),
	]
	subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def blast_query(query_fasta: Path, db_prefix: Path, max_evalue: float) -> list[BlastHit]:
	outfmt = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen"
	command = [
		"blastn",
		"-task",
		"blastn-short",
		"-query",
		str(query_fasta),
		"-db",
		str(db_prefix),
		"-dust",
		"no",
		"-soft_masking",
		"false",
		"-evalue",
		str(max_evalue),
		"-outfmt",
		outfmt,
	]
	result = subprocess.run(command, check=True, text=True, capture_output=True)
	hits: list[BlastHit] = []
	for line in result.stdout.splitlines():
		if not line.strip():
			continue
		fields = line.split("\t")
		hits.append(
			BlastHit(
				query_id=fields[0],
				subject_id=fields[1],
				pident=float(fields[2]),
				length=int(fields[3]),
				qstart=int(fields[4]),
				qend=int(fields[5]),
				sstart=int(fields[6]),
				send=int(fields[7]),
				evalue=float(fields[8]),
				bitscore=float(fields[9]),
				qlen=int(fields[10]),
				slen=int(fields[11]),
			)
		)
	return hits

def is_terminal_hit(hit: BlastHit, side: str, end_tolerance: int) -> bool:
	effective_tolerance = max(end_tolerance, 13)
	if side == "left":
		return hit.qstart <= effective_tolerance
	return hit.qend >= hit.qlen - effective_tolerance + 1

def select_best_hit(
	hits: Iterable[BlastHit],
	side: str,
	min_identity: float,
	min_query_coverage: float,
	min_alignment_length: int,
	end_tolerance: int,
) -> BlastHit | None:
	filtered = [
		hit
		for hit in hits
		if hit.pident >= min_identity
		and hit.qcov_pct >= min_query_coverage
		and hit.length >= min_alignment_length
		and is_terminal_hit(hit, side, end_tolerance)
	]
	if not filtered:
		return None
	filtered.sort(key=lambda hit: (hit.bitscore, hit.length, hit.pident), reverse=True)
	return filtered[0]

def classify_subject_coverage(hit: BlastHit, full_length: int = 180) -> str:
	tolerance = 3
	# Prioritize how much of the telomere start is missing, which is the critical part.
	missing_from_alignment_start = max(0, hit.subject_start - 1)
	missing_from_reference_start = max(0, full_length - hit.slen)
	missing_from_start = max(missing_from_alignment_start, missing_from_reference_start)

	if missing_from_start <= tolerance:
		return "Full"
	if missing_from_start <= 13 + tolerance:
		return "Truncated; missing palindrome I"
	if missing_from_start <= 40 + tolerance:
		return "Truncated; missing palindrome I and II"
	if missing_from_start <= 65 + tolerance:
		return "Truncated; missing palindrome I, II, and III"
	if missing_from_start <= 100 + tolerance:
		return "Truncated; missing palindrome I, II, and IV"
	return "No match"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description="Check assembly ends for Streptomycete telomere sequences using BLAST.",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	)
	
	parser.add_argument("assembly",
	type=Path,
	help="Input assembly FASTA")
	
	parser.add_argument("--telomeres",
	type=Path,
	help="Telomere sequence FASTA used to build the BLAST database",
	default=Path(__file__).parent / "telo_classes.fasta")
	
	parser.add_argument("--window-size",
	type=int,
	default=360,
	help="Terminal window size to scan from each contig end")
	
	parser.add_argument("--end-tolerance",
	"--terminal-bases",
	dest="end_tolerance",
	type=int,
	default=3,
	help="Allow a hit to start/end within this many bases of the contig end")
	
	parser.add_argument("--min-identity",
	type=float,
	default=70.0,
	help="Minimum percent identity for a hit")
	
	parser.add_argument("--min-query-coverage",
	type=float,
	default=0.0,
	help="Minimum percent query coverage for a hit")
	
	parser.add_argument("--min-alignment-length",
	type=int,
	default=50,
	help="Minimum aligned length for a hit")
	
	parser.add_argument("--max-evalue",
	type=float,
	default=1e-5,
	help="Maximum E-value to keep from BLAST")
	
	parser.add_argument("--full-table",
	action="store_true",
	help="Output full BLAST-detail table instead of compact output")
	
	parser.add_argument("--out",
	type=Path, default=None,
	help="Write results to this TSV file instead of stdout")
	return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
	args = parse_args(argv)

	if not args.assembly.exists():
		raise FileNotFoundError(f"Assembly FASTA not found: {args.assembly}")
	if not args.telomeres.exists():
		raise FileNotFoundError(f"Telomere FASTA not found: {args.telomeres}")

	records = list(read_fasta(args.assembly))
	if not records:
		raise ValueError(f"No sequences found in assembly FASTA: {args.assembly}")

	with tempfile.TemporaryDirectory(prefix="telotjek_") as temp_dir_name:
		temp_dir = Path(temp_dir_name)
		db_prefix = temp_dir / "telomeres_db"
		query_fasta = temp_dir / "assembly_ends.fasta"

		build_blast_db(args.telomeres, db_prefix)

		query_records: list[FastaRecord] = []
		query_meta: dict[str, tuple[str, str, int]] = {}
		for record in records:
			window = min(args.window_size, len(record.sequence))
			if window <= 0:
				continue
			left_id = f"{record.name}__left"
			right_id = f"{record.name}__right"
			query_records.append(FastaRecord(name=left_id, sequence=record.sequence[:window]))
			query_records.append(FastaRecord(name=right_id, sequence=record.sequence[-window:]))
			query_meta[left_id] = (record.name, "left", len(record.sequence))
			query_meta[right_id] = (record.name, "right", len(record.sequence))

		write_fasta(query_records, query_fasta)
		hits = blast_query(query_fasta, db_prefix, args.max_evalue)

		hits_by_query: dict[str, list[BlastHit]] = {}
		for hit in hits:
			hits_by_query.setdefault(hit.query_id, []).append(hit)

		output_handle = args.out.open("w", newline="") if args.out else sys.stdout
		close_output = args.out is not None
		try:
			writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
			if args.full_table:
				writer.writerow(
					[
						"Contig",
						"Side",
						"Class",
						"Coverage",
						"present",
						"query_length",
						"pident",
						"alignment_length",
						"subject_start",
						"subject_end",
						"query_coverage_pct",
						"query_start",
						"query_end",
						"strand",
						"evalue",
						"bitscore",
					]
				)
			else:
				writer.writerow(["Contig", "Side", "Class", "Coverage"])

			for query_id, (contig, side, query_length) in query_meta.items():
				best_hit = select_best_hit(
					hits_by_query.get(query_id, []),
					side=side,
					min_identity=args.min_identity,
					min_query_coverage=args.min_query_coverage,
					min_alignment_length=args.min_alignment_length,
					end_tolerance=args.end_tolerance,
				)
				if best_hit is None:
					if args.full_table:
						writer.writerow([contig, side, "", "No match", 0, query_length, "", "", "", "", "", "", "", "", "", ""])
					else:
						writer.writerow([contig, side, "", "No match"])
					continue

				coverage_class = classify_subject_coverage(best_hit)

				if args.full_table:
					writer.writerow(
						[
							contig,
							side,
							best_hit.subject_id,
							coverage_class,
							1,
							query_length,
							f"{best_hit.pident:.2f}",
							best_hit.length,
							best_hit.subject_start,
							best_hit.subject_end,
							f"{best_hit.qcov_pct:.2f}",
							best_hit.qstart,
							best_hit.qend,
							best_hit.strand,
							f"{best_hit.evalue:g}",
							f"{best_hit.bitscore:.1f}",
						]
					)
				else:
					writer.writerow([contig, side, best_hit.subject_id, coverage_class])
		finally:
			if close_output:
				output_handle.close()

	return 0


if __name__ == "__main__":
	raise SystemExit(main())
