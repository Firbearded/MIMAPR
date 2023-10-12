from scipy import interpolate
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def graph(k: str, t: int, h: int, width: int, heigth: int, data):
    data_gr = data[k][t][h]

    fig, ax = plt.subplots(figsize=(7, 7))
    im = ax.imshow(data_gr)
    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.5)

    fig.tight_layout()
    plt.show()


def graph_2(k: str, t: int, h: int, width: int, heigth: int, x, y, z):
    x1 = x[k][t][h]
    y1 = y[k][t][h]
    z1 = z[k][t][h]
    X, Y = np.meshgrid(x1, y1)
    f = interpolate.interp2d(x1, y1, z1)
    x2 = np.arange(h, width, 1)
    y2 = np.arange(h, heigth, 1)
    z2 = f(x2, y2)
    print(z2)

    fig, ax = plt.subplots(figsize=(7, 7))
    im = ax.imshow(z2)
    fig.tight_layout()
    plt.show()


def main():
    sim_t = 100
    width = 500
    heigth = 400
    h = [5, 10, 5, 10]
    data = {'exp': [], 'imp': []}
    x = {'exp': [], 'imp': []}
    y = {'exp': [], 'imp': []}
    z = {'exp': [], 'imp': []}
    files = ['explicit5.txt', 'explicit10.txt', 'implicit5.txt', 'implicit10.txt']
    file_k = ['exp', 'exp', 'imp', 'imp']
    for k in data.keys():
        for t in range(sim_t + 2):
            data[k].append({5: [], 10: []})
            x[k].append({5: [], 10: []})
            y[k].append({5: [], 10: []})
            z[k].append({5: [], 10: []})
            for h_ in h[:2]:
                for i in range(heigth // h_ + 1):
                    data[k][t][h_].append([])
                    for j in range(width // h_ + 1):
                        data[k][t][h_][i].append(-100)

    for f_i in range(len(files)):
        with open(files[f_i], 'r') as f:
            for line in f.readlines():
                s = line.split()
                x[file_k[f_i]][int(s[0])][h[f_i]].append(float(s[1]))
                y[file_k[f_i]][int(s[0])][h[f_i]].append(float(s[2]))
                z[file_k[f_i]][int(s[0])][h[f_i]].append(float(s[3]))
                data[file_k[f_i]][int(s[0])][h[f_i]][heigth // h[f_i] - int(s[2]) // h[f_i]][int(s[1]) // h[f_i]] = float(s[3])

    #print(data)
    graph('exp', 100, 10, width, heigth, data)
    graph('imp', 100, 10, width, heigth, data)
    graph('exp', 100, 5, width, heigth, data)
    graph('imp', 100, 5, width, heigth, data)
    #graph_2('exp', 100, 10, width, heigth, x, y, z)


if __name__ == "__main__":
    main()
