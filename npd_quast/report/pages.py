import os.path

from ..tools import SUPPORTED_TOOLS
import npd_quast.report.tables as tables
import npd_quast.report.metrics as metrics


EMPTY_REPORT = set()
RAW_REPORT = {'tool_answers.txt'}
FULL_REPORT = {
    'tool_answers.txt',
    'tool_page.html',
    'top_plot.png',
    'quantiles_plot.png',
}


class _AbstractPage:
    _s = None
    _npd_quast_folder = None
    _in = None

    def __init__(self, npd_quast_folder):
        self._npd_quast_folder = npd_quast_folder

    def _make_aside(self):
        aside = '''<aside>
<p class="p_aside">
    <a href="{0}">
        Total
    </a>
</p>'''.format(
            os.path.join('..', 'total_page.html') * self._in +
            'total_page.html' * (not self._in)
        )
        for report in os.listdir(
            os.path.join(
                self._npd_quast_folder, 'reports',
            )
        ):
            if os.path.isdir(
                os.path.join(
                    self._npd_quast_folder, 'reports', report,
                )
            ):
                if set(
                    os.listdir(
                        os.path.join(
                            self._npd_quast_folder, 'reports', report,
                        )
                    )
                ) != set():
                    aside += '''<p class="p_aside">
                        <a href="{0}">
                            {1}
                        </a>
                    </p>'''.format(
                        os.path.join('..', report, 'tool_page.html') * self._in +
                        os.path.join(report, 'tool_page.html') * (not self._in),
                        report,
                    )
        aside += '''<p class="p_aside">
    <a href="{0}">
        Help
    </a>
</p>
</aside>'''.format(
            os.path.join('..', 'about_metrics_page.html') * self._in +
            'about_metrics_page.html' * (not self._in)
        )
        return aside

    def __str__(self):
        return self._s


class TotalPage(_AbstractPage):
    _in = False

    def __init__(
            self,
            npd_quast_folder,
            true_answers,
            tool_answers_dict,
    ):
        super().__init__(npd_quast_folder)
        with open(
                os.path.join(
                    'npd_quast',
                    'report',
                    'templates',
                    'total_page.html',
                ),
        ) as total_page:
            s = total_page.read()
        self._s = s.replace(
            '$ASIDE$',
            super()._make_aside(),
        ).replace(
            '$TABLES$',
            ''.join(
                map(
                    str,
                    [
                        tables.TopTable(
                            true_answers,
                            tool_answers_dict,
                            tops=[1, 3, 5, 10],
                        ),
                        "<img src=\"top_plot.png\">",
                        tables.QuantilesTable(
                            true_answers,
                            tool_answers_dict,
                            quantiles=[25, 50, 75],
                        ),
                        "<img src=\"quantiles_plot.png\">",
                        tables.RankTable(
                            true_answers,
                            tool_answers_dict,
                        ),
                        tables.RRPTable(
                            true_answers,
                            tool_answers_dict,
                        ),
                        tables.MultiMedalsTable(
                            true_answers,
                            tool_answers_dict,
                        ),
                    ],
                ),
            ),
        )


class ToolPage(_AbstractPage):
    _in = True

    def __init__(
            self,
            npd_quast_folder,
            true_answers,
            tool_answers_dict,
            report,
    ):
        tool_answers = tool_answers_dict[report]
        super().__init__(npd_quast_folder)
        with open(
                os.path.join(
                    'npd_quast',
                    'report',
                    'templates',
                    'tool_page.html',
                ),
        ) as tool_page:
            s = tool_page.read()
        self._s = s.replace(
            '$TOOL_NAME$',
            report,
        ).replace(
            '$ASIDE$',
            super()._make_aside(),
        ).replace(
            '$MEAN_CORRECT_ANSWER_RANK$',
            str(metrics.mean_rank(true_answers, tool_answers)),
        ).replace(
            '$MEDIAN_CORRECT_ANSWER_RANK$',
            str(metrics.median_rank(true_answers, tool_answers)),
        ).replace(
            '$MEAN_RRP$',
            str(metrics.mean_rrp(true_answers, tool_answers)),
        ).replace(
            '$MEDIAN_RRP$',
            str(metrics.median_rrp(true_answers, tool_answers)),
        ).replace(
            '$MEAN_WEIGHTED_RRP$',
            str(metrics.mean_weighted_rrp(true_answers, tool_answers)),
        ).replace(
            '$MEDIAN_WEIGHTED_RRP$',
            str(metrics.median_weighted_rrp(true_answers, tool_answers)),
        )


class HelpPage(_AbstractPage):
    _in = False

    def __init__(self, npd_quast_folder):
        super().__init__(npd_quast_folder)
        with open(
                os.path.join(
                    'npd_quast',
                    'report',
                    'templates',
                    'help_page.html',
                )
        ) as about_metrics_page:
            s = about_metrics_page.read()
            self._s = s.replace(
                '$ASIDE$',
                super()._make_aside(),
            ).replace(
                '$SUPPORTED_TOOLS$',
                '<ul>' + '<br>'.join(
                    [
                        '<li>' + tool + '</li>'
                        for tool in SUPPORTED_TOOLS
                    ]
                ) + '</ul>'
            )
