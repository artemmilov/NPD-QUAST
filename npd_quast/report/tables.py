import npd_quast.report.metrics as metrics


class _AbstractTable:
    _title = None
    _metrics = None
    _single = False

    def __init__(self, true_answers, tool_answers_dict, single=False):
        self._true_answers = true_answers
        self._tool_answers_dict = tool_answers_dict
        self._single = single

    def __str__(self):
        table = '''
            <table border="1" align="center">
                    <caption>
                        <b>{0}</b>
                    </caption>
                    <tr>
        '''.format(self._title)
        for metric in [''] * (not self._single) + list(
                map(lambda m: m[0], self._metrics),
        ):
            table += '<th>{0}</th>\n'.format(metric)
        table += '</tr>\n'
        for tool, tool_answers in self._tool_answers_dict.items():
            if not self._single:
                table += '<tr>\n<th class="row_th">{0}</th>\n'.format(tool)
            for name, method in self._metrics:
                table += '<td>{0}</td>\n'.format(
                    method(
                        self._true_answers,
                        tool_answers,
                    ),
                )
            table += '</tr>\n'
        table += '</table>\n'
        return table


class TopTable(_AbstractTable):
    _title = 'Top'

    def __init__(self, true_answers, tool_answers_dict, tops, single=False):
        super().__init__(true_answers, tool_answers_dict, single)
        self._metrics = [
            (
                '{0}'.format(top),
                lambda tr_a, to_a:
                metrics.top_x(tr_a, to_a, top),
            )
            for top in tops
        ]


class QuantilesTable(_AbstractTable):
    _title = 'Quantiles'

    def __init__(self, true_answers, tool_answers_dict, quantiles, single=False):
        super().__init__(true_answers, tool_answers_dict, single)
        self._metrics = [
            (
                '{0}%'.format(quantile),
                lambda tr_a, to_a:
                metrics.k_quantile(tr_a, to_a, quantile),
            )
            for quantile in quantiles
        ]


class RankTable(_AbstractTable):
    _title = 'Correct answer rank'
    _metrics = [
        ('Mean', metrics.mean_rank),
        ('Median', metrics.median_rank),
    ]


class RRPTable(_AbstractTable):
    _title = 'RRP'
    _metrics = [
        ('Mean', metrics.mean_rrp),
        ('Median', metrics.median_rrp),
        ('Weighted mean', metrics.mean_weighted_rrp),
        ('Weighted median', metrics.median_weighted_rrp),
    ]


class MultiMedalsTable(_AbstractTable):
    _title = 'Medal score'
    _metrics = [
        ('Classic', metrics.classic_medal_score),
        ('F1', metrics.formula1_score),
        ('Gold', metrics.gold_medals),
        ('All', metrics.all_medals),
    ]

    def __str__(self):
        table = '''
            <table border="1" align="center">
                    <caption>
                        <b>{0}</b>
                    </caption>
                    <tr>
        '''.format(self._title)
        for metric in [''] + list(map(lambda m: m[0], self._metrics)):
            table += '<th>{0}</th>\n'.format(metric)
        table += '</tr>\n'
        for tool in self._tool_answers_dict:
            table += '<tr>\n<th class="row_th">{0}</th>\n'.format(tool)
            for name, method in self._metrics:
                table += '<td>{0}</td>\n'.format(
                    method(
                        self._true_answers,
                        self._tool_answers_dict,
                    )[tool],
                )
            table += '</tr>\n'
        table += '</table>\n'
        return table
