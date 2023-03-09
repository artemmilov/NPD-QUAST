import os.path
import shutil

from .pages import TotalPage, ToolPage, HelpPage
import npd_quast.report.plots as plots


def write_report(
        npd_quast_folder,
        tool_answers_dict,
        true_answers=None,
        true_answers_on=False,
        decoys_on=False
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
        if true_answers_on:
            print('1')
            print(tool_answers_dict)
            plots.write_interactive_top_plot(
                true_answers,
                {
                    report:
                        tool_answers_dict.get(report)
                },
                os.path.join(
                    npd_quast_folder,
                    'reports',
                    report,
                    'top_plot.html'
                ),
            )
            print('2')
            print(tool_answers_dict)
            plots.write_interactive_quantiles_plot(
                true_answers,
                {
                    report:
                    tool_answers_dict.get(report)
                },
                os.path.join(
                    npd_quast_folder,
                    'reports',
                    report,
                    'quantiles_plot.html'
                ),
            )
        # if decoys_on:
        #     plots.write_interactive_decoy_naive_method(
        #         tool_answers_dict[report],
        #         os.path.join(
        #             npd_quast_folder,
        #             'reports',
        #             report,
        #             'naive_method.html'
        #         ),
        #     )
        # with open(
        #     os.path.join(
        #         npd_quast_folder,
        #         'reports',
        #         report,
        #         'tool_page.html',
        #     ),
        #     'w',
        # ) as tool_page:
        #     tool_page.write(
        #         str(
        #             ToolPage(
        #                 npd_quast_folder,
        #                 true_answers,
        #                 tool_answers_dict,
        #                 report,
        #                 true_answers_on=true_answers_on,
        #                 decoys_on=decoys_on
        #             ),
        #         ),
        #     )
    # print(tool_answers_dict)
    with open(
        os.path.join(
            npd_quast_folder,
            'reports',
            'total_page.html',
        ),
        'w',
    ) as total_page:
        if true_answers_on:
            print('3')
            print(tool_answers_dict)
            plots.write_interactive_top_plot(
                true_answers,
                tool_answers_dict,
                os.path.join(
                    npd_quast_folder,
                    'reports',
                    'top_plot.html',
                ),
            )
            print('4')
            print(tool_answers_dict)
            plots.write_interactive_quantiles_plot(
                true_answers,
                tool_answers_dict,
                os.path.join(
                    npd_quast_folder,
                    'reports',
                    'quantiles_plot.html',
                ),
            )
    #     total_page.write(
    #         str(
    #             TotalPage(
    #                 npd_quast_folder,
    #                 true_answers,
    #                 tool_answers_dict,
    #                 true_answers_on=true_answers_on
    #             ),
    #         ),
    #     )
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
