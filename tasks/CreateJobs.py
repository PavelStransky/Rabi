with open("jobs.txt", "w") as f:
    for i in range(200):
        f.write(f"nice -20 julia111 Schnellbruder.jl {i} 7")
        f.write("\n")
    for i in range(200):
        f.write(f"nice -20 julia111 Schnellbruder.jl {i} 8")
        f.write("\n")
    for i in range(200):
        f.write(f"nice -20 julia111 Schnellbruder.jl {i} 9")
        f.write("\n")
