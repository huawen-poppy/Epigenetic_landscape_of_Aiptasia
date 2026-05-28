#!/usr/bin/env python

import argparse
from collections import defaultdict


def merge(read_prediction, site_prediction):
    """
    Convert TandemMod read-level prediction to site-level prediction.

    Expected input columns, no header:
    0 transcript_id
    1 transcriptome_site
    2 motif
    3 read_id
    4 label: mod / unmod
    5 probability
    """

    cutoffs = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]

    # key = transcript_id|site|motif
    # value = [mod_count_0.5, mod_count_0.6, ..., mod_count_0.95, total_reads]
    count_dict = defaultdict(lambda: [0, 0, 0, 0, 0, 0, 0])

    skipped = 0
    total_lines = 0

    with open(read_prediction, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            total_lines += 1
            items = line.split("\t")

            if len(items) < 6:
                skipped += 1
                continue

            transcript_id = items[0]
            site = items[1]
            motif = items[2]
            label = items[4]
            probability = float(items[5])

            site_id = f"{transcript_id}|{site}|{motif}"

            # total reads covering this site
            count_dict[site_id][6] += 1

            # count high-confidence modified reads only
            if label == "mod":
                for i, cutoff in enumerate(cutoffs):
                    if probability >= cutoff:
                        count_dict[site_id][i] += 1

    with open(site_prediction, "w") as out:
        out.write(
            "transcriptome_id\t"
            "site\t"
            "motif\t"
            "mod_reads_p0.5\t"
            "mod_reads_p0.6\t"
            "mod_reads_p0.7\t"
            "mod_reads_p0.8\t"
            "mod_reads_p0.9\t"
            "mod_reads_p0.95\t"
            "total_reads\t"
            "mod_rate_p0.5\t"
            "mod_rate_p0.6\t"
            "mod_rate_p0.7\t"
            "mod_rate_p0.8\t"
            "mod_rate_p0.9\t"
            "mod_rate_p0.95\n"
        )

        for site_id, counts in count_dict.items():
            transcript_id, site, motif = site_id.split("|")

            total_reads = counts[6]

            rates = [
                counts[i] / total_reads if total_reads > 0 else 0
                for i in range(6)
            ]

            out.write(
                f"{transcript_id}\t"
                f"{site}\t"
                f"{motif}\t"
                f"{counts[0]}\t"
                f"{counts[1]}\t"
                f"{counts[2]}\t"
                f"{counts[3]}\t"
                f"{counts[4]}\t"
                f"{counts[5]}\t"
                f"{total_reads}\t"
                f"{rates[0]:.6f}\t"
                f"{rates[1]:.6f}\t"
                f"{rates[2]:.6f}\t"
                f"{rates[3]:.6f}\t"
                f"{rates[4]:.6f}\t"
                f"{rates[5]:.6f}\n"
            )

    print(f"Finished.")
    print(f"Input lines processed: {total_lines}")
    print(f"Sites written: {len(count_dict)}")
    print(f"Skipped malformed lines: {skipped}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert TandemMod read-level predictions to site-level predictions."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input read-level TandemMod prediction file."
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output site-level prediction file."
    )

    args = parser.parse_args()

    merge(args.input, args.output)
