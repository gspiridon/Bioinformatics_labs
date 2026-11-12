from Bio import SeqIO
import matplotlib.pyplot as plt
from pathlib import Path

folder = Path(__file__).parent

genome1 = str(SeqIO.read(folder / "genome1.fasta", "fasta").seq).upper()
genome2 = str(SeqIO.read(folder / "genome2.fasta", "fasta").seq).upper()

print("Genome 1 length:", len(genome1), "bp")
print("Genome 2 length:", len(genome2), "bp")

enzymes = {
    "EcoRI": ("GAATTC", 1),
    "BamHI": ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1)
}

results = {}

for enzyme, (sequence, cut) in enzymes.items():
    print("\nEnzyme:", enzyme)
    print("Recognition sequence:", sequence)

    for name, dna in [("Genome 1", genome1), ("Genome 2", genome2)]:
        positions = []

        circular_dna = dna + dna[:len(sequence)-1]

        for i in range(len(dna)):
            if circular_dna[i:i+len(sequence)] == sequence:
                positions.append((i + cut) % len(dna))

        positions.sort()
        fragments = []

        if len(positions) == 0:
            fragments = [len(dna)]
        else:
            for i in range(len(positions)-1):
                fragments.append(positions[i+1] - positions[i])

            fragments.append(len(dna) - positions[-1] + positions[0])

        fragments.sort(reverse=True)
        results[(name, enzyme)] = fragments

        print(name)
        print("Restriction sites:", len(positions))
        print("Fragment sizes:", fragments)

    if results[("Genome 1", enzyme)] == results[("Genome 2", enzyme)]:
        print("Same restriction pattern")
    else:
        print("Different restriction pattern")

plt.figure(figsize=(11, 6))

lane = 1
labels = []

for enzyme in enzymes:
    for name in ["Genome 1", "Genome 2"]:
        fragments = results[(name, enzyme)]

        for size in fragments:
            plt.plot([lane - 0.3, lane + 0.3],
                     [size, size], color="blue", linewidth=2)

        labels.append(name + "\n" + enzyme)
        lane += 1

plt.xticks(range(1, 7), labels)
plt.ylabel("Fragment size (bp)")
plt.title("Gel electrophoresis - HPV16")
plt.xlim(0.5, 6.5)
plt.ylim(0, 8500)
plt.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.show()
