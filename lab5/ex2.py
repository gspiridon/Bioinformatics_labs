import itertools


def read_fasta(filename):
    sequence = ""

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if not line.startswith(">"):
                sequence += line

    return sequence.upper()


def generate_patterns():
    nucleotides = ["A", "C", "G", "T"]
    patterns = []

    for length in range(3, 7):
        for combination in itertools.product(
            nucleotides,
            repeat=length
        ):
            patterns.append("".join(combination))

    return patterns


def find_repetitions(sequence, patterns):
    results = []

    for pattern in patterns:
        for repeat_count in range(2, 4):
            repeated_sequence = pattern * repeat_count
            total_length = len(repeated_sequence)

            if 6 <= total_length <= 10:
                start = 0

                while True:
                    position = sequence.find(
                        repeated_sequence,
                        start
                    )

                    if position == -1:
                        break

                    results.append(
                        {
                            "motif": pattern,
                            "repetition": repeated_sequence,
                            "motif_length": len(pattern),
                            "total_length": total_length,
                            "start": position,
                            "end": position + total_length - 1
                        }
                    )

                    start = position + 1

    return results


def remove_overlapping_repetitions(results):
    results = sorted(
        results,
        key=lambda x: (
            x["start"],
            -x["total_length"]
        )
    )

    filtered_results = []

    for result in results:
        overlaps = False

        for saved_result in filtered_results:
            if not (
                result["end"] < saved_result["start"]
                or
                result["start"] > saved_result["end"]
            ):
                overlaps = True
                break

        if not overlaps:
            filtered_results.append(result)

    return filtered_results


def save_results(
    sequence,
    patterns,
    all_repetitions,
    filtered_repetitions
):
    with open(
        "dna_repetition_results.txt",
        "w"
    ) as file:

        file.write(
            f"DNA sequence length: "
            f"{len(sequence)} nucleotides\n"
        )

        file.write(
            f"Generated brute-force patterns: "
            f"{len(patterns)}\n"
        )

        file.write(
            f"All detected repetitions: "
            f"{len(all_repetitions)}\n"
        )

        file.write(
            f"Distinct non-overlapping repetitions: "
            f"{len(filtered_repetitions)}\n\n"
        )

        file.write(
            f"{'Motif':<10}"
            f"{'Repeated Sequence':<20}"
            f"{'Motif Length':<15}"
            f"{'Total Length':<15}"
            f"{'Start':<10}"
            f"{'End':<10}\n"
        )

        for result in filtered_repetitions:
            file.write(
                f"{result['motif']:<10}"
                f"{result['repetition']:<20}"
                f"{result['motif_length']:<15}"
                f"{result['total_length']:<15}"
                f"{result['start'] + 1:<10}"
                f"{result['end'] + 1:<10}\n"
            )


dna_sequence = read_fasta(
    "dna_sequence.fasta"
)

if not 1000 <= len(dna_sequence) <= 3000:
    raise ValueError(
        "The DNA sequence must contain between "
        "1000 and 3000 nucleotides."
    )


patterns = generate_patterns()

all_repetitions = find_repetitions(
    dna_sequence,
    patterns
)

filtered_repetitions = remove_overlapping_repetitions(
    all_repetitions
)


print(
    "DNA sequence length:",
    len(dna_sequence),
    "nucleotides"
)

print(
    "Generated brute-force patterns:",
    len(patterns)
)

print(
    "All detected repetitions:",
    len(all_repetitions)
)

print(
    "Distinct non-overlapping repetitions:",
    len(filtered_repetitions)
)

print()

print(
    f"{'Motif':<10}"
    f"{'Repeated Sequence':<20}"
    f"{'Motif Length':<15}"
    f"{'Total Length':<15}"
    f"{'Start':<10}"
    f"{'End':<10}"
)

for result in filtered_repetitions:
    print(
        f"{result['motif']:<10}"
        f"{result['repetition']:<20}"
        f"{result['motif_length']:<15}"
        f"{result['total_length']:<15}"
        f"{result['start'] + 1:<10}"
        f"{result['end'] + 1:<10}"
    )


save_results(
    dna_sequence,
    patterns,
    all_repetitions,
    filtered_repetitions
)

print()
print(
    "Results saved to "
    "dna_repetition_results.txt"
)