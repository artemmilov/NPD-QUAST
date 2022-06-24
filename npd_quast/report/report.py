import os.path
import shutil

from .pages import TotalPage, ToolPage, HelpPage
import npd_quast.report.plots as plots


def write_report(
        npd_quast_folder,
        true_answers,
        tool_answers_dict,
):
    for report in os.listdir(
        os.path.join(
            npd_quast_folder,
            'reports',
        )
    ):
        if not os.path.isdir(
            os.path.join(
                npd_quast_folder,
                'reports',
                report,
            )
        ):
            continue
        if set(
            os.listdir(
                os.path.join(
                    npd_quast_folder,
                    'reports',
                    report,
                )
            )
        ) == set():
            continue
        shutil.copy(
            os.path.abspath(
                os.path.join(
                    'npd_quast',
                    'report',
                    'templates',
                    'NPD-Quast.png',
                ),
            ),
            os.path.join(
                npd_quast_folder,
                'reports',
                'NPD-Quast.png',
            )
        )
        shutil.copy(
            os.path.abspath(
                os.path.join(
                    'npd_quast',
                    'report',
                    'templates',
                    'to_right.png',
                ),
            ),
            os.path.join(
                npd_quast_folder,
                'reports',
                'to_right.png',
            )
        )
        shutil.copy(
            os.path.abspath(
                os.path.join(
                    'npd_quast',
                    'report',
                    'templates',
                    'down.png',
                ),
            ),
            os.path.join(
                npd_quast_folder,
                'reports',
                'down.png',
            )
        )
        plots.write_top_plot(
            true_answers,
            {
                report:
                tool_answers_dict[report]
            },
            os.path.join(
                npd_quast_folder,
                'reports',
                report,
                'top_plot.png'
            ),
        )
        plots.write_quantiles_plot(
            true_answers,
            {
                report:
                tool_answers_dict[report]
            },
            os.path.join(
                npd_quast_folder,
                'reports',
                report,
                'quantiles_plot.png'
            ),
        )
        with open(
            os.path.join(
                npd_quast_folder,
                'reports',
                report,
                'tool_page.html',
            ),
            'w',
        ) as tool_page:
            tool_page.write(
                str(
                    ToolPage(
                        npd_quast_folder,
                        true_answers,
                        tool_answers_dict,
                        report,
                    ),
                ),
            )
    with open(
        os.path.join(
            npd_quast_folder,
            'reports',
            'total_page.html',
        ),
        'w',
    ) as total_page:
        plots.write_top_plot(
            true_answers,
            tool_answers_dict,
            os.path.join(
                npd_quast_folder,
                'reports',
                'top_plot.png',
            ),
        )
        plots.write_quantiles_plot(
            true_answers,
            tool_answers_dict,
            os.path.join(
                npd_quast_folder,
                'reports',
                'quantiles_plot.png',
            ),
        )
        total_page.write(
            str(
                TotalPage(
                    npd_quast_folder,
                    true_answers,
                    tool_answers_dict,
                ),
            ),
        )
    with open(
        os.path.join(
            npd_quast_folder,
            'reports',
            'about_metrics_page.html',
        ),
        'w',
    ) as about_metrics_page:
        about_metrics_page.write(
            str(HelpPage(npd_quast_folder)),
        )
