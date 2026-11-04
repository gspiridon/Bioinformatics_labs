import random

transposons = {
    "TE1": "ATGCGTACGTTAGC",
    "TE2": "TTCGAGCTAGCTAA",
    "TE3": "GGATCCGTTACGGA",
    "TE4": "CTAGGCTAACGTTC"
}

dna_length = 300
dna = list(''.join(random.choices("ATCG", k=dna_length)))

positions = []

for name, sequence in transposons.items():
    length = len(sequence)

    while True:
        start = random.randint(0, dna_length - length)
        end = start + length

        overlap = any(
            start < existing_end and end > existing_start
            for existing_start, existing_end in positions
        )

        if not overlap:
            positions.append((start, end))
            break

    dna[start:end] = list(sequence)

dna = ''.join(dna)

print("Artificial DNA sequence:")
print(dna)
print("\nDNA length:", len(dna), "bp")

print("\nDetected transposable elements:")

for name, sequence in transposons.items():
    start = dna.find(sequence)

    if start != -1:
        position = start + 1  
        length = len(sequence)
        end = position + length - 1

        print(
            f"{name}: Start = {position}, "
            f"End = {end}, Length = {length} bp"
        )
    else:
        print(f"{name}: Not detected")