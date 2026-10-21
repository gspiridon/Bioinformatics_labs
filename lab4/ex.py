import random
import matplotlib.pyplot as plt

FILENAME = "sequence.fasta"

genome = ""

with open(FILENAME, "r") as f:
    for line in f:
        if not line.startswith(">"):
            genome += line.strip()

samples = []

for _ in range(2000):
    length = random.randint(100, 200)
    start = random.randint(0, len(genome) - length)
    fragment = genome[start:start + length]
    samples.append((start, fragment))

samples.sort(key=lambda x: x[0])

start, first_fragment = samples[0]

reconstructed = first_fragment
current_end = start + len(first_fragment)

for start, fragment in samples[1:]:
    end = start + len(fragment)

    if start <= current_end:
        overlap = current_end - start

        if overlap < len(fragment):
            reconstructed += fragment[overlap:]
            current_end = end

    else:
        gap = start - current_end
        reconstructed += "N" * gap
        reconstructed += fragment
        current_end = end

coverage = [0] * len(genome)

for start, fragment in samples:
    for i in range(start, start + len(fragment)):
        coverage[i] += 1

uncovered = coverage.count(0)

coverage_percent = (
    (len(genome) - uncovered) / len(genome) * 100
)

average_depth = sum(coverage) / len(coverage)

lengths = [len(fragment) for _, fragment in samples]

zero_positions = [
    i for i, value in enumerate(coverage)
    if value == 0
]

print("Genome length:", len(genome))
print("Number of samples:", len(samples))
print("Minimum sample length:", min(lengths))
print("Maximum sample length:", max(lengths))
print("Average sample length:", round(sum(lengths) / len(lengths), 2))
print("Reconstructed sequence length:", len(reconstructed))
print("Uncovered bases:", uncovered)
print("Genome coverage:", round(coverage_percent, 2), "%")
print("Average coverage depth:", round(average_depth, 2), "x")
print("Uncovered positions:", zero_positions)

print("\nFirst 3 samples:")

for i in range(3):
    start, fragment = samples[i]

    print(
        f"Sample {i + 1}: position={start}, "
        f"length={len(fragment)}"
    )

    print(fragment[:60] + "...")

plt.figure(figsize=(12, 5))

plt.plot(range(len(coverage)), coverage)

plt.xlabel("Genome position")
plt.ylabel("Coverage depth")
plt.title("Genome Coverage")

plt.tight_layout()

plt.savefig("coverage_plot.png", dpi=300)

plt.show()

print("\nCoverage plot saved as coverage_plot.png")