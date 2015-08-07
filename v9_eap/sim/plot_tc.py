import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys
import re

def plot_platform(fname, ax, id):
    f = file(fname, 'r')
    start = []
    end = []
    for line in f.readlines():
        if "GW" in line:
            start = map(float, line.split()[3:6])
            end = map(float, line.split()[6:9])
            ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], color='black')
    f.close()

def main(fname):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    f = file(fname, 'r')
    lines = f.readlines()
    plt_f = re.findall("\"([^\"]+)\"", lines[0])
    plot_platform(plt_f[0], ax, id)
    colors = ['r','g','b','y','c','m','k','w','b','g']
    curr_color = -1
    x_p = []
    y_p = []
    z_p = []
    ants = 0
    for line in lines[1:]:
        if "add_antenna" in line:
            ants += 1
            curr_color += 1
            x_p.append([])
            y_p.append([])
            z_p.append([])
        elif "add_point" in line:
            p = re.findall("\(([^\"]+)\)", line)
            p = map(float, p[0].split(','))
            x_p[curr_color].append(p[0])
            y_p[curr_color].append(p[1])
            z_p[curr_color].append(p[2])
    assert len(x_p) == len(y_p) == len(z_p)
    tc_id = id
    for i in range(ants):
        ax.plot(x_p[i], y_p[i], z_p[i], "o", color = colors[i], label="Allowable placements for antenna " + str(i))
    f.close()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    #ax.set_title("Test Case %d" % id)
    labels = map(lambda x: "Antenna %d" % x, range(1,ants+1))
    leg = ax.legend(loc=0, numpoints=1)
    plt.savefig(fname.split('.')[0] + '.eps', format='eps', dpi=1000)
    plt.show()
    plt.clf()
main(sys.argv[1])
