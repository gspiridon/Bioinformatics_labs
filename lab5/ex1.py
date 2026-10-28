import random
import math
import matplotlib.pyplot as plt


def read_fasta(filename):
    sequence = ""

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if not line.startswith(">"):
                sequence += line

    return sequence.upper()


def molecular_weight(length_bp):
    return length_bp * 660


def migration_distance(length_bp):
    min_size = 100
    max_size = 3000

    max_distance = 90
    min_distance = 20

    distance = (
        max_distance
        +
        (
            math.log10(length_bp)
            -
            math.log10(min_size)
        )
        *
        (
            min_distance - max_distance
        )
        /
        (
            math.log10(max_size)
            -
            math.log10(min_size)
        )
    )

    return distance


def save_results(
    dna_length,
    lengths,
    weights,
    distances,
    sample_positions
):
    with open(
        "gel_electrophoresis_results.txt",
        "w"
    ) as file:

        file.write(
            f"Original DNA sequence length: "
            f"{dna_length} bp\n\n"
        )

        file.write("DNA fragments\n\n")

        for i in range(10):
            file.write(
                f"Fragment {i + 1}: "
                f"{lengths[i]} bp, "
                f"positions "
                f"{sample_positions[i][0] + 1}-"
                f"{sample_positions[i][1]}\n"
            )

        file.write("\nResults\n\n")

        file.write(
            f"{'Lane':<8}"
            f"{'Length (bp)':<15}"
            f"{'Molecular Weight (kDa)':<25}"
            f"{'Migration Distance (mm)':<25}\n"
        )

        for i in range(10):
            file.write(
                f"{i + 1:<8}"
                f"{lengths[i]:<15}"
                f"{weights[i] / 1000:<25.2f}"
                f"{distances[i]:<25.2f}\n"
            )


dna_sequence = read_fasta("dna_sequence.fasta")

print(
    "Original DNA sequence length:",
    len(dna_sequence),
    "bp"
)

if len(dna_sequence) < 1000:
    raise ValueError(
        "The DNA sequence must contain at least 1000 bases."
    )

max_fragment_size = min(
    3000,
    len(dna_sequence)
)

samples = []
sample_positions = []

for i in range(10):
    fragment_length = random.randint(
        100,
        max_fragment_size
    )

    start = random.randint(
        0,
        len(dna_sequence) - fragment_length
    )

    fragment = dna_sequence[
        start:start + fragment_length
    ]

    samples.append(fragment)

    sample_positions.append(
        (
            start,
            start + fragment_length
        )
    )


lengths = []
weights = []
distances = []

for fragment in samples:
    length = len(fragment)

    weight = molecular_weight(length)

    distance = migration_distance(length)

    lengths.append(length)
    weights.append(weight)
    distances.append(distance)


print()
print("DNA fragments")
print()

for i in range(10):
    print(
        f"Fragment {i + 1}: "
        f"{lengths[i]} bp, "
        f"positions "
        f"{sample_positions[i][0] + 1}-"
        f"{sample_positions[i][1]}"
    )


print()
print("Results")
print()

print(
    f"{'Lane':<8}"
    f"{'Length (bp)':<15}"
    f"{'Molecular Weight (kDa)':<25}"
    f"{'Migration Distance (mm)':<25}"
)

for i in range(10):
    print(
        f"{i + 1:<8}"
        f"{lengths[i]:<15}"
        f"{weights[i] / 1000:<25.2f}"
        f"{distances[i]:<25.2f}"
    )


save_results(
    len(dna_sequence),
    lengths,
    weights,
    distances,
    sample_positions
)


plt.figure(
    figsize=(10, 7)
)

for lane in range(1, 11):
    plt.plot(
        [lane - 0.3, lane + 0.3],
        [5, 5],
        color="black",
        linewidth=5
    )


for lane, distance, length in zip(
    range(1, 11),
    distances,
    lengths
):
    plt.plot(
        [lane - 0.3, lane + 0.3],
        [distance, distance],
        linewidth=6
    )

    plt.text(
        lane,
        distance + 3,
        f"{length} bp",
        ha="center",
        fontsize=9
    )


plt.xlim(
    0.5,
    10.5
)

plt.ylim(
    100,
    0
)

plt.xticks(
    range(1, 11),
    [
        "Lane 1",
        "Lane 2",
        "Lane 3",
        "Lane 4",
        "Lane 5",
        "Lane 6",
        "Lane 7",
        "Lane 8",
        "Lane 9",
        "Lane 10"
    ],
    rotation=45
)

plt.xlabel(
    "DNA Fragments"
)

plt.ylabel(
    "Migration Distance (mm)"
)

plt.title(
    "Gel Electrophoresis Simulation\n"
    "Short DNA Fragments Migrate Faster"
)

plt.grid(
    axis="y",
    alpha=0.3
)

plt.tight_layout()

plt.savefig(
    "gel_electrophoresis_output.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()