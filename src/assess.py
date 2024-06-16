import json
import matplotlib.pyplot as plt
import numpy as np

with open("mepr_slips.json", "r") as f:
    data = json.load(f)


x = np.arange(len(data))
scores = []
tx_rates = []

for i in data:
    scores.append(i["score"])
    tx_rates.append(i["tx_rate"])


fig, ax1 = plt.subplots()

# Plot y1 on the left y-axis
ax1.plot(x, scores, color='b', marker='o', label='y1')
ax1.set_xlabel('Position')
ax1.set_ylabel('Overlap scores', color='b')
ax1.tick_params(axis='y', labelcolor='b')

# Create a second y-axis
ax2 = ax1.twinx()

# Plot y2 on the right y-axis
ax2.plot(x, tx_rates, color='r', marker='o', label='y2')
ax2.set_ylabel('Tx rates', color='r')
ax2.tick_params(axis='y', labelcolor='r')

# Add a title
plt.title('Tx rates vs operator overlap scores for '+str(len(data))+' synthetic promoters')

# Show the plot
# fig.tight_layout()

print(data[22])

plt.show()