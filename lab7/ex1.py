
from Bio import SeqIO
from pathlib import Path

file_path = Path(__file__).parent / "sequence (1).fasta"

record = SeqIO.read(file_path, "fasta")
dna = str(record.seq).upper()

print("DNA length:", len(dna), "nucleotides")

enzymes = {
    "EcoRI": ("GAATTC", "G|AATTC"),
    "BamHI": ("GGATCC", "G|GATCC"),
    "HindIII": ("AAGCTT", "A|AGCTT")
}

for name, (recognition, cleavage) in enzymes.items():

    print("\nEnzyme:", name)
    print("Recognition sequence:", recognition)
    print("Cleavage site:", cleavage)

    positions = []

    for i in range(len(dna) - len(recognition) + 1):
        if dna[i:i + len(recognition)] == recognition:
            positions.append(i + 1)

    print("Number of restriction sites:", len(positions))

    if positions:
        print("Recognition site positions:", positions)
        print("Cut positions:", [pos + 1 for pos in positions])
    else:
        print("No restriction sites found")
