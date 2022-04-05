import npd_quast.metrics


def _str_best(true_answers, tool_answers):
    return ''.join(
        [
            'Correct in best {0}: {1:0.2f}%.\n'.format(
                best,
                npd_quast.metrics.percent_true_best(true_answers, tool_answers, best=best),
            )
            for best in [10, 100, 1000, None]
        ]
    )


def _str_medals(true_answers, tool_answers):
    return ''.join(
        [
            'Classic medal score: {0}.\n'.format(
                npd_quast.metrics.classic_medal_score(true_answers, tool_answers),
            ),
            'F1 score: {0}.\n'.format(
                npd_quast.metrics.formula1_score(true_answers, tool_answers),
            ),
            'Gold medals: {0}.\n'.format(
                npd_quast.metrics.gold_medals(true_answers, tool_answers),
            ),
            'All medals: {0}.\n'.format(
                npd_quast.metrics.all_medals(true_answers, tool_answers),
            ),
        ]
    )


def _str_ranks(true_answers, tool_answers, default_ranks=None):
    if default_ranks is None:
        default_ranks = []
    return ''.join(
        [
            ''.join(
                [
                    'Mean rank: {0:0.2f}. Default rank is {1}.\n'.format(
                        npd_quast.metrics.mean_rank(true_answers, tool_answers, default_rank),
                        default_rank,
                    ),
                    'Median rank: {0:0.2f}. Default rank is {1}.\n'.format(
                        npd_quast.metrics.median_rank(true_answers, tool_answers, default_rank),
                        default_rank,
                    ),
                ],
            )
            for default_rank in [None] + default_ranks
        ],
    )


def _str_rprs(true_answers, tool_answers):
    return ''.join(
        [
            'Mean RPR: {0:0.2f}.\n'.format(
                npd_quast.metrics.mean_rpr(true_answers, tool_answers),
            ),
            'Median RPR: {0:0.2f}.\n'.format(
                npd_quast.metrics.median_rpr(true_answers, tool_answers),
            ),
            'Mean weighted RPR: {0:0.2f}.\n'.format(
                npd_quast.metrics.mean_weighted_rpr(true_answers, tool_answers),
            ),
            'Median weighted RPR: {0:0.2f}.\n'.format(
                npd_quast.metrics.median_weighted_rpr(true_answers, tool_answers),
            ),
        ],
    )


def _str_quantiles(true_answers, tool_answers):
    return ''.join(
        [
            '{0}% quantile: {1:0.2f}.\n'.format(
                k,
                npd_quast.metrics.k_quantile(true_answers, tool_answers, k),
            )
            for k in [25, 50, 75, 100]
        ],
    )


def write_report(report_folder, true_answers, tool_answers, default_ranks=None):
    with open(report_folder, 'w') as report:
        report.write(
            ''.join(
                [
                    _str_best(true_answers, tool_answers),
                    _str_medals(true_answers, tool_answers),
                    _str_ranks(true_answers, tool_answers, default_ranks),
                    _str_rprs(true_answers, tool_answers),
                    _str_quantiles(true_answers, tool_answers),
                ],
            ),
        )


def _str_multi_medals(true_answers, tool_answers_dict):
    return '\n'.join(
        [
            'Classic multi medal score:\n {0}'.format(
                ''.join(
                    '\t{0}:\t{1}\n'.format(tool_title, score)
                    for tool_title, score in
                    npd_quast.metrics.classic_multi_medal_score(
                        true_answers,
                        tool_answers_dict,
                    ).items()
                ),
            ),
            'F1 multi score:\n {0}'.format(
                ''.join(
                    '\t{0}:\t{1}\n'.format(tool_title, score)
                    for tool_title, score in
                    npd_quast.metrics.formula1_multi_score(
                        true_answers,
                        tool_answers_dict,
                    ).items()
                ),
            ),
            'Gold multi medals:\n {0}'.format(
                ''.join(
                    '\t{0}:\t{1}\n'.format(tool_title, score)
                    for tool_title, score in
                    npd_quast.metrics.gold_multi_medals(
                        true_answers,
                        tool_answers_dict,
                    ).items()
                ),
            ),
            'All multi medals:\n {0}'.format(
                ''.join(
                    '\t{0}:\t{1}\n'.format(tool_title, score)
                    for tool_title, score in
                    npd_quast.metrics.all_multi_medals(
                        true_answers,
                        tool_answers_dict,
                    ).items()
                ),
            ),
        ]
    )


def write_multi_report(multi_report_folder, true_answers, tool_answers_dict):
    with open(multi_report_folder, 'w') as report:
        report.write(
            ''.join(
                [
                    _str_multi_medals(true_answers, tool_answers_dict),
                ],
            ),
        )
