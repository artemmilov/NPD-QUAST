from matplotlib.ticker import MultipleLocator

from .metrics import top_x, k_quantile, fdr
import matplotlib.pyplot as plt

import plotly.graph_objs as go

import numpy as np

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


def write_interactive_top_plot(true_answers, tool_answers_dict, folder):
    # объявляем фигуру
    fig = go.Figure()

    n = 10
    #print(tool_answers_dict)
    for i, tool in enumerate(tool_answers_dict.keys()):
        tool_answers = tool_answers_dict[tool]
        if sum(map(len, tool_answers.values())) > n:
            n = sum(map(len, tool_answers.values()))
        #print(tool_answers)
        tops = [
            top_x(true_answers, tool_answers, top)
            for top in range(1, n)
        ]

        # добавляем график
        fig.add_trace(
            go.Scatter(
                x=np.arange(1, n), y=np.array(tops),  # данные
                name=tool,  # имя в легенде
                marker=dict(color='#00CC66'),  # цвет в html-формате
                opacity=0.8,  # прозрачность
                line={'width': 3}  # свойства линии - толщина
            )
        )
        # if max(tops) > m:
        #     m = max(tops)
        # if len(tool_answers_dict) == 1:
        #     ax.plot(range(1, n), tops, alpha=1.0)
        # else:
        #     ax.plot(range(1, n), tops, alpha=1.0, color=COLORS[i], label=tool)
        #     legend = True

    # свойства фигуры
    fig.update_layout(
        height=450, width=700,  # размер фигуры
        title_text='Top x',  # заголовок графика
        title_font_size=16,  # размер заголовка
        plot_bgcolor='rgba(0,0,0,0.05)',  # цвет фона
    )

    # параметры оси абсцисс
    fig.update_xaxes(
        range=[-1.5, 1.5],  # ограничение графика
        zeroline=True,  # рисовать линию x=0
        zerolinewidth=2  # толщина линии x=0
    )

    # параметры оси ординат
    fig.update_yaxes(
        zeroline=True,  # рисовать линию y=0
        zerolinewidth=2,  # толщина линии y=0
        zerolinecolor='LightGray'  # цвет линии y=0
    )

    # показать график
    fig.write_html(folder)


def write_interactive_quantiles_plot(true_answers, tool_answers_dict, folder):
    print('a')
    print(tool_answers_dict)
    # print(tool_answers_dict)
    # объявляем фигуру
    fig = go.Figure()

    n = 10
    for i, tool in enumerate(tool_answers_dict.keys()):
        tool_answers = tool_answers_dict[tool]
        quantiles = [
            k_quantile(true_answers, tool_answers, k)
            for k in range(10, 90)
        ]
        # добавляем график
        fig.add_trace(
            go.Scatter(
                x=np.arange(10, 90), y=np.array(quantiles),  # данные
                name=tool,  # имя в легенде
                marker=dict(color='#00CC66'),  # цвет в html-формате
                opacity=0.8,  # прозрачность
                line={'width': 3}  # свойства линии - толщина
            )
        )

    # свойства фигуры
    fig.update_layout(
        height=450, width=700,  # размер фигуры
        title_text='Quantiles',  # заголовок графика
        title_font_size=16,  # размер заголовка
        plot_bgcolor='rgba(0,0,0,0.05)',  # цвет фона
    )

    # параметры оси абсцисс
    fig.update_xaxes(
        range=[-1.5, 1.5],  # ограничение графика
        zeroline=True,  # рисовать линию x=0
        zerolinewidth=2  # толщина линии x=0
    )

    # параметры оси ординат
    fig.update_yaxes(
        zeroline=True,  # рисовать линию y=0
        zerolinewidth=2,  # толщина линии y=0
        zerolinecolor='LightGray'  # цвет линии y=0
    )

    # показать график
    fig.write_html(folder)

    print('b')
    print(tool_answers_dict)


def write_interactive_decoy_naive_method(tool_answers, folder):
    # mass_spectra = []
    # for challenge in os.listdir(os.path.join(folder, 'challenges')):
    #     for specter in os.listdir(os.path.join(folder, 'challenges', challenge, 'spectra')):
    #         f = os.path.join(folder, 'challenges', challenge, 'spectra', specter)
    #         mass_spectra.append(MassSpecter(f))


    # объявляем фигуру
    fig = go.Figure()

    # Question! Tool_answers
    n = sum([len(x) for x in tool_answers.values() if x != []])
    #print(n)
    fdrs = [
        fdr(tool_answers, k)
        for k in range(1, n + 1)
    ]

    # добавляем график
    fig.add_trace(
        go.Scatter(
            x=np.arange(1, n), y=np.array(fdrs),  # данные
            name='Naive decoys',  # имя в легенде
            marker=dict(color='#00CC66'),  # цвет в html-формате
            opacity=0.8,  # прозрачность
            line={'width': 3}  # свойства линии - толщина
        )
    )

    # свойства фигуры
    fig.update_layout(
        height=450, width=700,  # размер фигуры
        title_text='Naive decoys',  # заголовок графика
        title_font_size=16,  # размер заголовка
        plot_bgcolor='rgba(0,0,0,0.05)',  # цвет фона
    )

    # параметры оси абсцисс
    fig.update_xaxes(
        range=[-1.5, 1.5],  # ограничение графика
        zeroline=True,  # рисовать линию x=0
        zerolinewidth=2  # толщина линии x=0
    )

    # параметры оси ординат
    fig.update_yaxes(
        zeroline=True,  # рисовать линию y=0
        zerolinewidth=2,  # толщина линии y=0
        zerolinecolor='LightGray'  # цвет линии y=0
    )

    # показать график
    fig.write_html(folder)
