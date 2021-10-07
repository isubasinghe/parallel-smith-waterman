import seaborn as sns
import matplotlib.pyplot as plt

x = [8, 8, 16, 16, 24, 24, 32, 32]
y = [1000, 1100, 400, 500, 350, 250, 105, 125]

x = [1,2,8,16,24,32,36]
y = [923129129, 491227097, 169576852, 119517599, 112753964, 111237046, 108483498]

eff = [1451965419/m for m in y]
for (i, p) in enumerate(x):
    eff[i] = eff[i]

print(eff)

sns.set_theme(color_codes=True)
def main():
    assert(len(x) == len(y))
    sns_fig = sns.scatterplot(x=x, y=eff)
    fig = sns_fig.get_figure()
    plt.xlabel("Threads")
    plt.ylabel("Efficiency")
    fig.savefig("eff.png")

if __name__ == '__main__':
    main()
