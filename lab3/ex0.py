from collections import Counter
import matplotlib.pyplot as plt
import subprocess
import re


GENETIC_CODE = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "Stop", "TAG": "Stop",
    "TGT": "Cys", "TGC": "Cys", "TGA": "Stop", "TGG": "Trp",

    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",

    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",

    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
}


def read_fasta(filename):
    sequence = ""

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if not line.startswith(">"):
                sequence += line.upper()

    return sequence


def reverse_complement(sequence):
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }

    return "".join(
        complement.get(base, "N")
        for base in reversed(sequence)
    )


def count_codons_in_frame(sequence, frame):
    counts = Counter()

    for i in range(frame, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]

        if all(base in "ATCG" for base in codon):
            counts[codon] += 1

    return counts


def count_all_frames(sequence):
    total_counts = Counter()

    for frame in range(3):
        total_counts.update(
            count_codons_in_frame(sequence, frame)
        )

    reverse = reverse_complement(sequence)

    for frame in range(3):
        total_counts.update(
            count_codons_in_frame(reverse, frame)
        )

    return total_counts


def show_top_codons(counts, number=10):
    top = counts.most_common(number)

    print("\nTop codons:")

    for codon, frequency in top:
        amino_acid = GENETIC_CODE[codon]

        print(
            f"{codon:4} "
            f"{amino_acid:5} "
            f"{frequency}"
        )

    return top


def plot_codons(counts, number=10):
    top = counts.most_common(number)

    codons = [codon for codon, frequency in top]
    frequencies = [frequency for codon, frequency in top]

    plt.figure(figsize=(10, 6))
    plt.bar(codons, frequencies)

    plt.xlabel("Codon")
    plt.ylabel("Frequency")
    plt.title(f"Top {number} most frequent codons")

    plt.tight_layout()
    plt.savefig("top_codons.png", dpi=300)
    plt.close()

    print("\nChart saved as top_codons.png")


def create_llm_prompt(top_codons):
    amino_acids = []

    for codon, frequency in top_codons:
        amino_acid = GENETIC_CODE[codon]

        if amino_acid != "Stop" and amino_acid not in amino_acids:
            amino_acids.append(amino_acid)

    amino_acid_text = ", ".join(amino_acids)

    prompt = f"""
The amino acids associated with the 10 most frequent codons are:

{amino_acid_text}

Create exactly one markdown table with common prepared dishes
that contain relatively low amounts of each amino acid.

Use exactly these columns:

Amino Acid | Dish 1 | Dish 2 | Dish 3

Rules:
- Include each amino acid exactly once.
- Include exactly 3 prepared dishes for each amino acid.
- Use complete dish names, not single ingredients.
- Do not use vague items such as "fruit", "vegetables", "bread",
  "soup", "cereal", "chicken breast", "rice", or "tofu" by themselves.
- Each entry should describe a recognizable meal or prepared dish.
- Good examples of the expected level of detail are:
  "vegetable soup with potatoes and carrots",
  "pasta with tomato and basil sauce",
  "boiled rice with steamed vegetables",
  "baked potato with vegetable filling",
  "oatmeal with apple and cinnamon",
  "vegetable risotto",
  "polenta with roasted vegetables",
  "mashed potatoes with steamed carrots",
  "vegetable pasta salad",
  "rice porridge with fruit".
- Choose dishes that are relatively low in the specified amino acid.
- Do not exclude any food category automatically.
- You may include meat, fish, dairy, grains, vegetables, fruits,
  or mixed dishes if they are appropriate.
- Prefer common dishes that people actually eat.
- Use varied dishes across the table.
- Do not repeat the same dish too often.
- Do not invent unusual or artificial food varieties.
- Do not include stop codons.
- Do not list codons as amino acids.
- Do not add notes or explanations outside the table.

Return only the table.
"""

    return prompt


def clean_ansi(text):
    ansi_escape = re.compile(
        r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])'
    )

    return ansi_escape.sub('', text)


def ask_ollama(prompt, model="llama3.2"):
    ollama_path = r"C:\Users\Alina\AppData\Local\Programs\Ollama\ollama.exe"

    try:
        result = subprocess.run(
            [ollama_path, "run", model],
            input=prompt,
            text=True,
            capture_output=True,
            encoding="utf-8",
            errors="replace",
            timeout=120
        )

        if result.returncode != 0:
            print("\nOllama returned an error:")
            print(result.stderr)
            return None

        return result.stdout

    except FileNotFoundError:
        print("\nOllama executable was not found.")
        return None

    except subprocess.TimeoutExpired:
        print("\nOllama took too long to respond.")
        return None


def main():
    filename = "GCA_008369605.1_ASM836960v1_genomic.fna"

    print("Reading genome...")

    genome = read_fasta(filename)

    print(f"Genome length: {len(genome):,} bp")

    print("Counting codons in all 6 reading frames...")

    codon_counts = count_all_frames(genome)

    sorted_codons = codon_counts.most_common()

    print("\nAll codons sorted by frequency:")

    for codon, frequency in sorted_codons:
        print(
            codon,
            GENETIC_CODE[codon],
            frequency
        )

    top10 = show_top_codons(codon_counts, 10)

    plot_codons(codon_counts, 10)

    prompt = create_llm_prompt(top10)

    print("\nGenerated LLM prompt:")
    print(prompt)

    print("\nSending prompt to Ollama...")

    answer = ask_ollama(prompt)

    if answer:
        clean_answer = clean_ansi(answer)

        print("\nLLM response:")
        print(clean_answer)

        with open("llm_foods.txt", "w", encoding="utf-8") as file:
            file.write(clean_answer)

        print("\nLLM response saved as llm_foods.txt")


if __name__ == "__main__":
    main()