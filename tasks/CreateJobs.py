with open("jobs.txt", "w") as f:
    for i in range(90):
        f.write(f"nice -20 julia111 Schnellbruder.jl {i} 0")
        f.write("\n")
    for i in range(90):
        f.write(f"nice -20 julia111 Schnellbruder.jl {i} 2")
        f.write("\n")
