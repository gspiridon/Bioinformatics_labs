import math
import matplotlib.pyplot as plt


def calculate_salt_adjusted_tm(sequence: str, na_conc: float = 0.05) -> dict:
    """Calculates the melting temperature (Tm) of a DNA sequence using the salt-adjusted formula:

    Tm = 81.5 + 16.6 * log10([Na+]) + 0.41 * (%GC) - 600 / length

    Parameters:
    - sequence (str): DNA sequence string.
    - na_conc (float): Sodium ion concentration [Na+] in Molar
      (default: 0.05 M = 50 mM).
    """

    seq = sequence.upper().strip()

    valid_bases = {"A", "T", "G", "C"}
    invalid_chars = set(seq) - valid_bases

    if invalid_chars:
        raise ValueError(
            f"Invalid nucleotide(s) found: {', '.join(invalid_chars)}"
        )

    length = len(seq)

    if length == 0:
        raise ValueError("Sequence length must be greater than 0.")

    if na_conc <= 0:
        raise ValueError("Sodium concentration [Na+] must be greater than 0.")

    counts = {base: seq.count(base) for base in ["A", "T", "G", "C"]}

    gc_count = counts["G"] + counts["C"]
    percent_gc = (gc_count / length) * 100

    tm = (
        81.5
        + (16.6 * math.log10(na_conc))
        + (0.41 * percent_gc)
        - (600 / length)
    )

    print(f"Sequence: {seq}")
    print(f"Length: {length} bp")
    print(f"Base Counts: {counts}")
    print(f"%GC Content: {percent_gc:.2f}%")
    print(f"[Na+] Concentration: {na_conc} M")
    print(f"Melting Temperature (Tm): {tm:.2f}°C\n")

    plt.figure(figsize=(7, 5))

    colors = ["#2ca02c", "#d62728", "#1f77b4", "#ff7f0e"]
    bars = plt.bar(counts.keys(), counts.values(), color=colors)

    plt.xlabel("Nucleotide Base", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)

    plt.title(
        f'Nucleotide Frequencies & Salt-Adjusted Tm ({tm:.2f}°C)\n'
        f'Sequence: "{seq}" | [Na+] = {na_conc} M',
        fontsize=11,
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

    chart_filename = "salt_adjusted_tm_chart.png"
    plt.savefig(chart_filename, dpi=300, bbox_inches="tight")

    print(f"Chart saved as: {chart_filename}")

    plt.show()

    return {
        "sequence": seq,
        "length": length,
        "counts": counts,
        "percent_gc": percent_gc,
        "tm_celsius": tm,
        "chart_file": chart_filename,
    }


if __name__ == "__main__":
    sequence = "TAACGCTA"
    calculate_salt_adjusted_tm(sequence, na_conc=0.05)