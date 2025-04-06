with open("jobs.txt", "w") as f:
    for i in range(90):
        f.write(f"nice 30 julia18 Schnellbruder.jl 0 {i}")
        f.write("\n")
    for i in range(90):
        f.write(f"nice 30 julia18 Schnellbruder.jl 2 {i}")
        f.write("\n")
