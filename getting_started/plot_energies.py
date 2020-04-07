import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


def get_energies(file):
    '''Get energies per iteration from an output file.'''
    energies = []
    with open(file,'r') as file:
        for line in file:
            tmp = line.rstrip().split(' E= ')
            if len(tmp) >= 2:
                tmp = tmp[1].split('  ')[0]
                energies.append(float(tmp))
    return energies


def main():
    e_sto3g = get_energies('h2o_sto3g.out')
    e_dzp = get_energies('h2o_dzp.out')
    fig, ax = plt.subplots()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.plot(e_sto3g, label='STO-3G')
    plt.plot(e_dzp, label='DZP')
    plt.xlabel('iteration')
    plt.ylabel('energies [Ha]')
    plt.legend(loc='center right')
    plt.show()


if __name__ == "__main__":
    main()
