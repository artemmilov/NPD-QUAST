import metrics
from collections import defaultdict


def main():
    answers_folder = input()
    true_answers = {}
    tool_answers = defaultdict(list)
    with open(answers_folder + '/true_answers.txt') as true_answers_data:
        for true_answer in true_answers_data.read().split('\n'):
            if true_answer != '':
                true_answers[int(true_answer.split('\t')[0])] = true_answer.split('\t')[1]
    with open(answers_folder + '/tool_answers.txt') as tool_answers_data:
        for tool_answer in tool_answers_data.read().split('\n'):
            if tool_answer != '':
                tool_answers[int(tool_answer.split('\t')[0])].append(
                    (
                        tool_answer.split('\t')[1],
                        tool_answer.split('\t')[2],
                    ),
                )
    for best in [10, 100, 1000, None]:
        print(
            'Correct in best {0}: {1:0.2f}%.'.format(
                best,
                metrics.percent_true_best(true_answers, tool_answers, best=best),
            ),
        )
    print(
        'Classic medal score: {0}.'.format(
            metrics.classic_medal_score(true_answers, tool_answers),
        ),
    )
    print(
        'F1 score: {0}.'.format(
            metrics.formula1_score(true_answers, tool_answers),
        ),
    )
    print(
        'Gold medals: {0}.'.format(
            metrics.gold_medals(true_answers, tool_answers),
        ),
    )
    print(
        'All medals: {0}.'.format(
            metrics.all_medals(true_answers, tool_answers),
        ),
    )
    print(
        'Mean rank: {0:0.2f}.'.format(
            metrics.mean_rank(true_answers, tool_answers),
        ),
    )
    print(
        'Median rank: {0:0.2f}.'.format(
            metrics.median_rank(true_answers, tool_answers),
        ),
    )
    print(
        'Mean RPR: {0:0.2f}.'.format(
            metrics.mean_rpr(true_answers, tool_answers),
        ),
    )
    print(
        'Median RPR: {0:0.2f}.'.format(
            metrics.median_rpr(true_answers, tool_answers),
        ),
    )
    print(
        'Mean weighted RPR: {0:0.2f}.'.format(
            metrics.mean_weighted_rpr(true_answers, tool_answers),
        ),
    )
    print(
        'Median weighted RPR: {0:0.2f}.'.format(
            metrics.median_weighted_rpr(true_answers, tool_answers),
        ),
    )
    for k in [25, 50, 75, 100]:
        print(
            '{0}% quantile: {1:0.2f}.'.format(
                k,
                metrics.k_quantile(true_answers, tool_answers, k),
            )
        )


if __name__ == '__main__':
    main()  # converter_from_NPDTools/test_main_input
