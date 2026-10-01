import matplotlib.pyplot as plt
cells = 100
hours = [0]
population = [cells]

for hour in range(1, 11):
    cells = cells * 2
    hours.append(hour)
    population.append(cells)

plt.plot(hours, population)
plt.xlabel("Time (hours)")
plt.ylabel("Cells")
plt.show()
