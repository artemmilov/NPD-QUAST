from numpy import mean, median, quantile


def percent_true_best(true_answers, tool_answers, best=None):
    sorted_answers = sorted(
        [
            [scan, tool_answer[0], tool_answer[1]]
            for scan, scan_answers in tool_answers.items()
            for tool_answer in scan_answers
        ],
        key=lambda ans: ans[2],
    )
    total = len(sorted_answers)
    if best is None:
        best = total
    correct_matches = 0
    for answer in sorted_answers[:best]:
        scan, tool_inchi, _ = answer
        correct_matches += (tool_inchi == true_answers[scan])
    return correct_matches / min(best, total) * 100


def _abstract_medal_score(true_answers, tool_answers, medals):
    score = 0
    for scan, true_inchi in true_answers.items():
        sorted_tool_answers = sorted(
            tool_answers[scan],
            key=lambda ans: ans[1],
        )
        sorted_inches = tuple(
            map(
                lambda ans: ans[0],
                sorted_tool_answers,
            )
        )
        if true_inchi in sorted_inches:
            pos = sorted_inches.index(true_inchi)
            if pos < len(medals):
                score += medals[pos]
    return score


def classic_medal_score(true_answers, tool_answers):
    return _abstract_medal_score(true_answers, tool_answers, [5, 3, 1])


def formula1_score(true_answers, tool_answers):
    return _abstract_medal_score(
        true_answers,
        tool_answers,
        [24, 18, 15, 12, 10, 8, 6, 4, 2, 1],
    )


def gold_medals(true_answers, tool_answers):
    return _abstract_medal_score(true_answers, tool_answers, [1])


def all_medals(true_answers, tool_answers):
    return _abstract_medal_score(true_answers, tool_answers, [1, 1, 1])


def _get_sorted_inches_for_scan(tool_answers, scan):
    sorted_answers = sorted(
        [
            [tool_answer[0], tool_answer[1]]
            for tool_answer in tool_answers[scan]
        ],
        key=lambda ans: ans[1],
    )
    result = []
    last_scan = None
    for sorted_answer in sorted_answers:
        if sorted_answer[1] == last_scan:
            result[-1].append(sorted_answer[0])
        else:
            result.append([sorted_answer[0]])
            last_scan = sorted_answer[1]
    return result


def _get_ranks(true_answers, tool_answers, default_rank=None):
    ranks = []
    for scan, true_inchi in true_answers.items():
        sorted_inches = _get_sorted_inches_for_scan(tool_answers, scan)
        found = False
        for i, set_inches in enumerate(sorted_inches):
            if true_inchi in set_inches:
                before = sum(map(len, sorted_inches[:i]))
                ans_rank = mean(
                    list(
                        map(
                            lambda x: before + x,
                            range(0, len(set_inches)),
                        ),
                    ),
                )
                ranks.append(ans_rank)
                found = True
        if not found:
            if default_rank is not None:
                ranks.append(default_rank)
    return ranks


def mean_rank(true_answers, tool_answers, default_rank=None):
    return mean(_get_ranks(true_answers, tool_answers, default_rank))


def median_rank(true_answers, tool_answers, default_rank=None):
    return median(_get_ranks(true_answers, tool_answers, default_rank))


def _get_rprs(true_answers, tool_answers):
    rprs = []
    for scan, true_inchi in true_answers.items():
        sorted_inches = _get_sorted_inches_for_scan(tool_answers, scan)
        for i, set_inches in enumerate(sorted_inches):
            if true_inchi in set_inches:
                before = sum(map(len, sorted_inches[:i]))
                after = sum(map(len, sorted_inches[i + 1:]))
                total = sum(list(map(len, sorted_inches)))
                if total > 1:
                    rprs.append(1 / 2 * (1 - (before - after) / (total - 1)))
                else:
                    rprs.append(1 / 2)
    return rprs


def mean_rpr(true_answers, tool_answers):
    return mean(_get_rprs(true_answers, tool_answers))


def median_rpr(true_answers, tool_answers):
    return median(_get_rprs(true_answers, tool_answers))


def _get_weighted_rprs(true_answers, tool_answers):
    weighted_rprs = []
    for scan, true_inchi in true_answers.items():
        sorted_answers = sorted(
            [
                [tool_answer[0], tool_answer[1]]
                for tool_answer in tool_answers[scan]
            ],
            key=lambda ans: ans[1],
            reverse=True,
        )
        sorted_inches = list(
            map(
                lambda ans: ans[0],
                sorted_answers,
            ),
        )
        sorted_scores = list(
            map(
                lambda ans: float(ans[1]),
                sorted_answers,
            ),
        )
        if true_inchi in sorted_inches:
            true_pos = sorted_inches.index(true_inchi)
            true_score = sorted_scores[true_pos]
            up_normalized = sum(score for score in sorted_scores if score > true_score) / sum(sorted_scores)
            same_normalized = sum(score for score in sorted_scores if score == true_score) / sum(sorted_scores)
            weighted_rprs.append(1 - up_normalized - same_normalized)
    return weighted_rprs


def mean_weighted_rpr(true_answers, tool_answers):
    return mean(_get_weighted_rprs(true_answers, tool_answers))


def median_weighted_rpr(true_answers, tool_answers):
    return median(_get_weighted_rprs(true_answers, tool_answers))


def k_quantile(true_answers, tool_answers, k=50):
    positions = []
    for scan, true_inchi in true_answers.items():
        sorted_inches = _get_sorted_inches_for_scan(tool_answers, scan)
        for i, set_inches in enumerate(sorted_inches):
            if true_inchi in set_inches:
                before = sum(map(len, sorted_inches[:i]))
                ans_rank = mean(
                    list(
                        map(
                            lambda x: before + x,
                            range(0, len(set_inches)),
                        ),
                    ),
                )
                positions.append(ans_rank)
    if len(positions) != 0:
        return quantile(positions, k / 100)
    return 0
