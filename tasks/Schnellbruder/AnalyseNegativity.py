import numpy as np
import matplotlib.pyplot as plt

# Read numbers from file
# Replace 'data.txt' with your filename

PATH = 'd:/results/Rabi/schnellbruder/'

J2 = 4
OMEGA = 1
R = 20

LAMBDA = -np.sqrt(2) / 5
DELTA = 0.5

negativity = []
purity = []

times = np.linspace(0, 400, 4000, endpoint=False)

for time in times:
    try:
        data = np.loadtxt(f'{PATH}Wigner_Rabi(2j={J2}, ω={OMEGA}, R={R}, λ={LAMBDA}, δ={DELTA})_{time:.2f}.txt')

    except Exception as e:
        print(f"Error loading {time:.2f}: {e}")
        continue

    y = data[:,2]
    w = sum(y)
    y /= w

    negativity.append(0.5 * (sum(abs(y)) - 1))
    purity.append(sum(y * y))

    print(time)

try:
    datax = np.loadtxt(f'{PATH}Jx_Rabi(2j={J2}, ω={OMEGA}, R={R}, λ={LAMBDA}, δ={DELTA}).txt')[:,1]
    datay = np.loadtxt(f'{PATH}Jy_Rabi(2j={J2}, ω={OMEGA}, R={R}, λ={LAMBDA}, δ={DELTA}).txt')[:,1]
    dataz = np.loadtxt(f'{PATH}Jz_Rabi(2j={J2}, ω={OMEGA}, R={R}, λ={LAMBDA}, δ={DELTA}).txt')[:,1]

    ts = np.loadtxt(f'{PATH}Jx_Rabi(2j={J2}, ω={OMEGA}, R={R}, λ={LAMBDA}, δ={DELTA}).txt')[:,0]

except Exception as e:
    print(f"Error loading file: {e}")
        
data = np.sqrt(datax * datax + datay * datay + dataz * dataz) / (J2 / 2)

plt.figure(figsize=(15, 10))
plt.plot(times, negativity, label='Negativity')
plt.plot(ts, data, label="Purity")
plt.ylim(0, 1)
plt.show()

np.savetxt(f'{PATH}Negativity_Rabi(2j={J2}, ω={OMEGA}, R={R}, λ={LAMBDA}, δ={DELTA}).txt', np.column_stack((times, negativity)))
np.savetxt(f'{PATH}Purity_Rabi(2j={J2}, ω={OMEGA}, R={R}, λ={LAMBDA}, δ={DELTA}).txt', np.column_stack((ts, data)))