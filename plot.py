import seaborn as sns
import matplotlib.pyplot as plt

x = [8, 8, 16, 16, 24, 24, 32, 32]
y = [1000, 1100, 400, 500, 350, 250, 105, 125]
sns.set_theme(color_codes=True)
def main():
    assert(len(x) == len(y))
    sns_fig = sns.regplot(x=x, y=y, ci=None)
    fig = sns_fig.get_figure()
    plt.xlabel("Threads")
    plt.ylabel("Microseconds")
    fig.savefig("res.png")

if __name__ == '__main__':
    main()
