from matplotlib.ticker import MultipleLocator

from .metrics import top_x, k_quantile
import matplotlib.pyplot as plt

COLORS = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']


def write_top_plot(true_answers, tool_answers_dict, folder):
    legend = False
    n = 10
    m = 10
    fig, ax = plt.subplots()
    ax.set_title('Top x')
    ax.set_xlabel('x')
    ax.set_ylabel('')
    for i, tool in enumerate(tool_answers_dict.keys()):
        tool_answers = tool_answers_dict[tool]
        if sum(map(len, tool_answers.values())) > n:
            n = sum(map(len, tool_answers.values()))
        tops = [
            top_x(true_answers, tool_answers, top)
            for top in range(1, n)
        ]
        if max(tops) > m:
            m = max(tops)
        if len(tool_answers_dict) == 1:
            ax.plot(range(1, n), tops, alpha=1.0)
        else:
            ax.plot(range(1, n), tops, alpha=1.0, color=COLORS[i], label=tool)
            legend = True

    if n > 5:
        ax.xaxis.set_major_locator(MultipleLocator((n // 5 + 9) // 10 * 10))
        ax.xaxis.set_minor_locator(MultipleLocator((n // 5 + 9) // 10))
    else:
        ax.xaxis.set_major_locator(MultipleLocator(2))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
    if m > 10:
        ax.yaxis.set_major_locator(MultipleLocator((m // 10 + 9) // 10 * 10))
        ax.yaxis.set_minor_locator(MultipleLocator((m // 10 + 9) // 10))
    else:
        ax.yaxis.set_major_locator(MultipleLocator(2))
        ax.yaxis.set_minor_locator(MultipleLocator(1))

    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    ax.grid(True)
    if legend:
        ax.legend()

    fig.savefig(folder)


def write_quantiles_plot(true_answers, tool_answers_dict, folder):
    legend = False
    m = 10
    fig, ax = plt.subplots()
    ax.set_title('Quantile k%')
    ax.set_xlabel('k')
    ax.set_ylabel('')
    for i, tool in enumerate(tool_answers_dict.keys()):
        tool_answers = tool_answers_dict[tool]
        quantiles = [
            k_quantile(true_answers, tool_answers, k)
            for k in range(10, 90)
        ]
        if max(quantiles) - min(quantiles) > m:
            m = max(quantiles) - min(quantiles)
        if len(tool_answers_dict) == 1:
            ax.plot(range(10, 90), quantiles, alpha=1.0)
        else:
            ax.plot(range(10, 90), quantiles, alpha=1.0, color=COLORS[i], label=tool)
            legend = True

    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(2))
    if m > 1:
        ax.yaxis.set_major_locator(MultipleLocator(m / 5))
        ax.yaxis.set_minor_locator(MultipleLocator(m / 50))
    else:
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    ax.grid(True)
    if legend:
        ax.legend()

    fig.savefig(folder)
