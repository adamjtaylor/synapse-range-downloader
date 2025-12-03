#!/usr/bin/env python3
"""
FASTQ Composition Detection using Sourmash

Detects organism composition in FASTQ files using k-mer sketching and
containment metrics. Useful for contamination detection and species verification.

Usage:
    python fastq_composition.py sample.fastq.gz
    python fastq_composition.py --url https://example.com/file.fastq.gz
    python fastq_composition.py sample.fastq.gz --ref-dir ./references
"""

import sourmash
import requests
import gzip
import os
import json
import sys
import argparse
from pathlib import Path


class CompositionAgent:
    """Agent for detecting organism composition in FASTQ files."""

    def __init__(self, ref_dir="references/", ksize=31, scaled=1000):
        """
        Initialize the composition agent.

        Args:
            ref_dir: Directory containing reference .sig files
            ksize: K-mer size for sketching (default: 31)
            scaled: Scaling factor for MinHash (default: 1000)
        """
        self.ref_dir = Path(ref_dir)
        self.ksize = ksize
        self.scaled = scaled
        self.refs = {}

        if not self.ref_dir.exists():
            raise FileNotFoundError(
                f"Reference directory '{ref_dir}' not found. "
                f"Create it with: mkdir {ref_dir}"
            )

        self._load_references()

    def _load_references(self):
        """Load all reference signatures from the reference directory."""
        sig_files = list(self.ref_dir.glob("*.sig"))

        if not sig_files:
            raise FileNotFoundError(
                f"No .sig files found in {self.ref_dir}. "
                f"Please add reference signatures."
            )

        print(f"Loading references from {self.ref_dir}...", file=sys.stderr)

        for sig_path in sig_files:
            try:
                # Load signature with matching ksize
                sig = sourmash.load_one_signature(sig_path, ksize=self.ksize)

                # Extract name from filename
                name = sig_path.stem.replace('_', ' ').title()
                self.refs[name] = sig

                print(f"  Loaded: {name}", file=sys.stderr)
            except Exception as e:
                print(f"  Warning: Could not load {sig_path.name}: {e}",
                      file=sys.stderr)

        if not self.refs:
            raise ValueError("No valid reference signatures could be loaded")

        print(f"Agent ready with {len(self.refs)} reference(s)\n",
              file=sys.stderr)

    def inspect_file(self, fastq_path, read_limit=50000):
        """
        Analyze composition of a local FASTQ file.

        Args:
            fastq_path: Path to FASTQ file (can be gzipped)
            read_limit: Maximum number of reads to process

        Returns:
            dict: Composition analysis results
        """
        fastq_path = Path(fastq_path)

        if not fastq_path.exists():
            return {"error": f"File not found: {fastq_path}"}

        print(f"Analyzing: {fastq_path.name}", file=sys.stderr)
        print(f"Sampling up to {read_limit:,} reads...\n", file=sys.stderr)

        # Create sketch
        sketch = sourmash.MinHash(n=0, ksize=self.ksize, scaled=self.scaled)

        # Process FASTQ file
        reads_seen = 0
        open_func = gzip.open if fastq_path.suffix == '.gz' else open

        try:
            with open_func(fastq_path, 'rt') as f:
                batch = []
                for line in f:
                    batch.append(line.strip())
                    if len(batch) == 4:
                        # Extract sequence (line 2 of FASTQ record)
                        seq = batch[1]
                        sketch.add_sequence(seq, force=True)
                        batch = []
                        reads_seen += 1

                        if reads_seen >= read_limit:
                            break
        except Exception as e:
            return {"error": f"Failed to process file: {str(e)}"}

        if reads_seen == 0:
            return {"error": "No reads found in file"}

        # Calculate composition
        return self._calculate_composition(
            sketch=sketch,
            reads_sampled=reads_seen,
            source=str(fastq_path)
        )

    def inspect_url(self, url, read_limit=50000):
        """
        Analyze composition of a FASTQ file from a URL.

        Args:
            url: URL to FASTQ file (must be gzipped)
            read_limit: Maximum number of reads to process

        Returns:
            dict: Composition analysis results
        """
        print(f"Streaming from URL...", file=sys.stderr)
        print(f"Sampling up to {read_limit:,} reads...\n", file=sys.stderr)

        # Create sketch
        sketch = sourmash.MinHash(n=0, ksize=self.ksize, scaled=self.scaled)

        # Stream from URL
        reads_seen = 0

        try:
            with requests.get(url, stream=True) as r:
                r.raise_for_status()

                with gzip.open(r.raw, 'rt') as f:
                    batch = []
                    for line in f:
                        batch.append(line.strip())
                        if len(batch) == 4:
                            seq = batch[1]
                            sketch.add_sequence(seq, force=True)
                            batch = []
                            reads_seen += 1

                            if reads_seen >= read_limit:
                                break
        except Exception as e:
            return {"error": f"Failed to stream from URL: {str(e)}"}

        if reads_seen == 0:
            return {"error": "No reads found in stream"}

        # Calculate composition
        return self._calculate_composition(
            sketch=sketch,
            reads_sampled=reads_seen,
            source=url
        )

    def _calculate_composition(self, sketch, reads_sampled, source):
        """
        Calculate organism composition by comparing sketch to references.

        Args:
            sketch: MinHash sketch of the sample
            reads_sampled: Number of reads processed
            source: Source identifier (file path or URL)

        Returns:
            dict: Composition analysis with contamination flags
        """
        if len(sketch) == 0:
            return {
                "error": "Sketch is empty - no valid k-mers found",
                "source": source
            }

        composition = {}
        total_explained = 0.0

        # Calculate containment for each reference
        for name, ref_sig in self.refs.items():
            # Containment: fraction of sample k-mers found in reference
            common = sketch.count_common(ref_sig.minhash)
            containment = common / len(sketch) if len(sketch) > 0 else 0

            # Filter noise (below 1% is likely sequencing error or homology)
            if containment > 0.01:
                composition[name] = round(containment, 4)
                total_explained += containment

        # Calculate unknown fraction
        unknown = max(0, 1.0 - total_explained)

        # Determine contamination status
        is_mixed = len(composition) > 1
        phix_contamination = composition.get('Phix', 0) > 0.05
        cross_contamination = is_mixed and not phix_contamination
        is_contaminated = phix_contamination or cross_contamination

        # Generate warning message
        warning = None
        if phix_contamination:
            warning = "PhiX control contamination detected - possible sequencing failure"
        elif cross_contamination:
            species_list = "/".join(composition.keys())
            warning = f"{species_list} mixture detected"
        elif unknown > 0.5:
            warning = "High unknown content - organism may not be in reference set"

        # Build result
        result = {
            "source": source,
            "reads_sampled": reads_sampled,
            "composition_estimate": composition,
            "unknown_content": round(unknown, 4),
            "is_mixed": is_mixed,
            "is_contaminated": is_contaminated,
        }

        if warning:
            result["contamination_warning"] = warning

        return result


def main():
    """Main entry point for the composition detection tool."""
    parser = argparse.ArgumentParser(
        description="Detect organism composition in FASTQ files using sourmash",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze a local file
  python fastq_composition.py sample.fastq.gz

  # Stream from URL
  python fastq_composition.py --url https://example.com/file.fastq.gz

  # Use custom reference directory
  python fastq_composition.py sample.fastq.gz --ref-dir ./my_refs

  # Limit reads processed
  python fastq_composition.py sample.fastq.gz --reads 10000
        """
    )

    parser.add_argument(
        'fastq',
        nargs='?',
        help='Path to FASTQ file (can be gzipped)'
    )

    parser.add_argument(
        '--url',
        help='URL to FASTQ file (must be gzipped)'
    )

    parser.add_argument(
        '--ref-dir',
        default='references/',
        help='Directory containing reference .sig files (default: references/)'
    )

    parser.add_argument(
        '--reads',
        type=int,
        default=50000,
        help='Maximum number of reads to process (default: 50000)'
    )

    parser.add_argument(
        '--ksize',
        type=int,
        default=31,
        help='K-mer size (default: 31)'
    )

    parser.add_argument(
        '--scaled',
        type=int,
        default=1000,
        help='Scaling factor for MinHash (default: 1000)'
    )

    args = parser.parse_args()

    # Validate input
    if not args.fastq and not args.url:
        parser.error("Either FASTQ file path or --url must be provided")

    if args.fastq and args.url:
        parser.error("Cannot specify both file path and --url")

    try:
        # Initialize agent
        agent = CompositionAgent(
            ref_dir=args.ref_dir,
            ksize=args.ksize,
            scaled=args.scaled
        )

        # Analyze composition
        if args.url:
            result = agent.inspect_url(args.url, read_limit=args.reads)
        else:
            result = agent.inspect_file(args.fastq, read_limit=args.reads)

        # Output JSON
        print(json.dumps(result, indent=2))

        # Exit with error code if there was a problem
        if "error" in result:
            sys.exit(1)

        # Exit with error code if contamination detected
        if result.get("is_contaminated"):
            sys.exit(2)

    except FileNotFoundError as e:
        print(json.dumps({"error": str(e)}, indent=2))
        sys.exit(1)

    except Exception as e:
        print(json.dumps({"error": f"Unexpected error: {str(e)}"}, indent=2),
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
