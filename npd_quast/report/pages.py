import os.path

from ..tools import SUPPORTED_TOOLS
import npd_quast.report.tables as tables


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
        for tool in SUPPORTED_TOOLS:
            tool_name = SUPPORTED_TOOLS[tool]().name()
            if tool_name in os.listdir(
                os.path.join(
                    self._npd_quast_folder, 'reports'
                )
            ):
                aside += '''<p class="p_aside">
    <a href="{0}">
        {1}
    </a>
</p>'''.format(
                    os.path.join('..', tool_name, 'tool_page.html') * self._in +
                    os.path.join(tool_name, 'tool_page.html') * (not self._in),
                    tool,
                )
        aside += '''<aside>
<p class="p_aside">
    <a href="{0}">
        About metrics
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
                        tables.SingleMedalsTable(
                            true_answers,
                            tool_answers_dict,
                        ),
                        tables.QuantilesTable(
                            true_answers,
                            tool_answers_dict,
                            quantiles=[25, 50, 70, 100],
                        ),
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
            tool,
    ):
        tool_answers = {tool: tool_answers_dict[tool]}
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
            tool,
        ).replace(
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
                            tool_answers,
                            tops=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                        ),
                        tables.SingleMedalsTable(
                            true_answers,
                            tool_answers,
                        ),
                        tables.QuantilesTable(
                            true_answers,
                            tool_answers,
                            quantiles=[5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100],
                        ),
                        tables.RankTable(
                            true_answers,
                            tool_answers,
                        ),
                        tables.RRPTable(
                            true_answers,
                            tool_answers,
                        ),
                    ],
                ),
            ),
        )


class AboutMetricsPage(_AbstractPage):
    _in = False

    def __init__(self, npd_quast_folder):
        super().__init__(npd_quast_folder)
        with open(
                os.path.join(
                    'npd_quast',
                    'report',
                    'templates',
                    'about_metrics_page.html',
                )
        ) as about_metrics_page:
            s = about_metrics_page.read()
            self._s = s.replace(
                '$ASIDE$',
                super()._make_aside(),
            )
