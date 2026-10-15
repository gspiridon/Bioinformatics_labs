from collections import Counter
import matplotlib.pyplot as plt

S = "ATCGACGTAACGTA"

pairs = []

for i in range(len(S) - 1):
    pairs.append(S[i:i+2])

freq = Counter(pairs)

print(freq)

plt.bar(freq.keys(), freq.values())
plt.xlabel("combinations")
plt.ylabel("frequency")
plt.title("DNA frequencies")
plt.show()