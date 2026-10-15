import matplotlib.pyplot as plt


def analyze_dna_sequence(sequence: str) -> dict:
    """Calculates the melting temperature (Tm) of a DNA sequence using the
    Wallace rule (Tm = 4*(G+C) + 2*(A+T)) and plots the nucleotide frequencies.
    """

    seq = sequence.upper().strip()

    valid_bases = {"A", "T", "G", "C"}
    invalid_chars = set(seq) - valid_bases

    if invalid_chars:
        raise ValueError(
            f"Invalid nucleotide(s) found: {', '.join(invalid_chars)}"
        )

    if len(seq) == 0:
        raise ValueError("Sequence length must be greater than 0.")

    counts = {base: seq.count(base) for base in ["A", "T", "G", "C"]}

    gc_count = counts["G"] + counts["C"]
    at_count = counts["A"] + counts["T"]

    tm = 4 * gc_count + 2 * at_count

    print(f"Sequence: {seq}")
    print(f"Length: {len(seq)} bp")
    print(f"Base Counts: {counts}")
    print(f"Melting Temperature (Tm): {tm}°C")

    # Generate Frequency Chart
    plt.figure(figsize=(7, 5))

    colors = ["#2ca02c", "#d62728", "#1f77b4", "#ff7f0e"]
    bars = plt.bar(counts.keys(), counts.values(), color=colors)

    plt.xlabel("Nucleotide Base", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)

    plt.title(
        f'Nucleotide Frequencies & Melting Temp (Tm = {tm}°C)\n'
        f'Sequence: "{seq}"',
        fontsize=12,
    )

    plt.ylim(0, max(counts.values()) + 1)

    for bar in bars:
        height = bar.get_height()

        plt.text(
            bar.get_x() + bar.get_width() / 2,
            height + 0.1,
            f"{int(height)}",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=11,
        )

    plt.tight_layout()

    chart_filename = "wallace_tm_chart.png"
    plt.savefig(chart_filename, dpi=300, bbox_inches="tight")

    print(f"Chart saved as: {chart_filename}")

    plt.show()

    return {
        "sequence": seq,
        "length": len(seq),
        "counts": counts,
        "tm_celsius": tm,
        "chart_file": chart_filename,
    }


if __name__ == "__main__":
    dna_seq = "TAACGCTA"
    analyze_dna_sequence(dna_seq)