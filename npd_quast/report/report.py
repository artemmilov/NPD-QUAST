import os.path

from .pages import TotalPage, ToolPage, AboutMetricsPage
from ..tools import SUPPORTED_TOOLS


def write_report(
        npd_quast_folder,
        true_answers,
        tool_answers_dict,
):
    for tool_name, tool in SUPPORTED_TOOLS.items():
        if tool().name() not in os.listdir(
            os.path.join(
                npd_quast_folder,
                'reports',
            )
        ):
            continue
        with open(
            os.path.join(
                npd_quast_folder,
                'reports',
                tool().name(),
                'tool_page.html'.format(tool().name()),
            ),
            'w',
        ) as tool_page:
            tool_page.write(
                str(
                    ToolPage(
                        npd_quast_folder,
                        true_answers,
                        tool_answers_dict,
                        tool_name,
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
            str(AboutMetricsPage(npd_quast_folder)),
        )
