with open("jobs.txt", "w") as f:
    for i in range(200):
        f.write(f"nice -20 julia111 Schnellbruder.jl {i} 10")
        f.write("\n")
    for i in range(200):
        f.write(f"nice -20 julia111 Schnellbruder.jl {i} 11")
        f.write("\n")
    for i in range(200):
        f.write(f"nice -20 julia111 Schnellbruder.jl {i} 12")
        f.write("\n")
