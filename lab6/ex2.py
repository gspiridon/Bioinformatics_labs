with open("sequence.fasta", "r") as file:
    dna = ""

    for line in file:
        if not line.startswith(">"):
            dna += line.strip().upper()

print("DNA length:", len(dna), "bp")

max_length = 10
min_length = 4

found = set()

print("\nPossible transposable elements:")

for length in range(max_length, min_length - 1, -1):
    for i in range(len(dna) - length + 1):

        sequence = dna[i:i + length]

        if sequence in found:
            continue

        positions = []

        for j in range(len(dna) - length + 1):
            if dna[j:j + length] == sequence:
                positions.append(j + 1)

        if len(positions) >= 2:
            found.add(sequence)

            for start in positions:
                end = start + length - 1

                print(
                    f"Sequence: {sequence}, "
                    f"Start: {start}, "
                    f"End: {end}, "
                    f"Length: {length} bp"
                )

if not found:
    print("No repeated sequences found.")