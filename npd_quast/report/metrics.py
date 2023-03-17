from numpy import mean, median, quantile

ROUND = 2


def top_x(true_answers, tool_answers, x=None):

    all_answers = []
    for challenge, challenge_answers in tool_answers.items():
        for challenge_answer in challenge_answers:
            all_answers.append(
                (
                    challenge,
                    challenge_answer[0],
                    challenge_answer[1],
                )
            )
    all_answers = sorted(
        all_answers,
        key=lambda ans: ans[2],
    )
    total = len(all_answers)
    if (x is None) or (x > total):
        x = total
    all_answers = all_answers[:x]
    return len(list(filter(
        lambda ans: ans[1] == true_answers[ans[0]],
        all_answers,
    )))


def _get_sorted_inches_for_scan(tool_answers, scan):
    result = []
    if tool_answers.get(scan) is None:
        return result
    sorted_answers = sorted(
        [
            [tool_answer[0], tool_answer[1]]
            for tool_answer in tool_answers.get(scan)       # !!!!!!!!!!!!!!!!!!!!!!!
        ],
        key=lambda ans: ans[1],
    )
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
                try:
                    ans_rank = mean(
                        list(
                            map(
                                lambda x: before + x,
                                range(0, len(set_inches)),
                            ),
                        ),
                    )
                except RuntimeWarning:
                    print('_get_ranks')
                ranks.append(ans_rank)
                found = True
        if not found:
            if default_rank is not None:
                ranks.append(default_rank)
    return ranks


def mean_rank(true_answers, tool_answers, default_rank=None):
    try:
        return round(
            float(mean(_get_ranks(true_answers, tool_answers, default_rank))),
            ROUND,
        )
    except RuntimeWarning:
        print('mean_rank')


def median_rank(true_answers, tool_answers, default_rank=None):
    return round(
        float(
            median(_get_ranks(true_answers, tool_answers, default_rank))
        ),
        ROUND,
    )


def _get_rrps(true_answers, tool_answers):
    rrps = []
    for scan, true_inchi in true_answers.items():
        sorted_inches = _get_sorted_inches_for_scan(tool_answers, scan)
        for i, set_inches in enumerate(sorted_inches):
            if true_inchi in set_inches:
                before = sum(map(len, sorted_inches[:i]))
                after = sum(map(len, sorted_inches[i + 1:]))
                total = sum(list(map(len, sorted_inches)))
                if total > 1:
                    rrps.append(1 / 2 * (1 - (before - after) / (total - 1)))
                else:
                    rrps.append(1 / 2)
    return rrps


def mean_rrp(true_answers, tool_answers):
    try:
        return round(
            float(mean(_get_rrps(true_answers, tool_answers))),
            ROUND,
        )
    except RuntimeWarning:
        print('mean_rrp')


def median_rrp(true_answers, tool_answers):
    try:
        return round(
            float(median(_get_rrps(true_answers, tool_answers))),
            ROUND,
        )
    except RuntimeWarning:
        print('median_rrp')


def _get_weighted_rrps(true_answers, tool_answers):
    weighted_rrps = []
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
            weighted_rrps.append(1 - up_normalized - same_normalized)
    return weighted_rrps


def mean_weighted_rrp(true_answers, tool_answers):
    try:
        return round(
            float(
                mean(_get_weighted_rrps(true_answers, tool_answers))
            ),
            ROUND,
        )
    except RuntimeWarning:
        print('mean_weighted_rrp')


def median_weighted_rrp(true_answers, tool_answers):
    return round(
        float(
            median(_get_weighted_rrps(true_answers, tool_answers))
        ),
        ROUND,
    )


def k_quantile(true_answers, tool_answers, k=50):
    positions = []
    for scan, true_inchi in true_answers.items():
        sorted_inches = _get_sorted_inches_for_scan(tool_answers, scan)
        for i, set_inches in enumerate(sorted_inches):
            if true_inchi in set_inches:
                before = sum(map(len, sorted_inches[:i]))
                try:
                    ans_rank = mean(
                        list(
                            map(
                                lambda x: before + x,
                                range(0, len(set_inches)),
                            ),
                        ),
                    )
                except RuntimeWarning:
                    print('k_quantile')
                positions.append(ans_rank)
    if len(positions) != 0:
        return round(quantile(positions, k / 100), ROUND)
    return 0


def _abstract_medal_score(true_answers, tool_answers_dict, medals):
    scores_dict = {tool_title: 0 for tool_title in tool_answers_dict}
    for scan, true_inchi in true_answers.items():
        positions_dict = {}
        for tool_title, tool_answers in tool_answers_dict.items():
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
                positions_dict[tool_title] = pos
        for place in range(0, min(len(medals), len(positions_dict))):
            tool_title = sorted(
                positions_dict.items(),
                key=lambda result: result[1],
            )[place][0]
            scores_dict[tool_title] += medals[place]
    return scores_dict


def classic_medal_score(true_answers, tool_answers_dict):
    return _abstract_medal_score(true_answers, tool_answers_dict, [5, 3, 1])


def formula1_score(true_answers, tool_answers_dict):
    return _abstract_medal_score(
        true_answers,
        tool_answers_dict,
        [24, 18, 15, 12, 10, 8, 6, 4, 2, 1],
    )


def gold_medals(true_answers, tool_answers_dict):
    return _abstract_medal_score(true_answers, tool_answers_dict, [1])


def all_medals(true_answers, tool_answers_dict):
    return _abstract_medal_score(true_answers, tool_answers_dict, [1, 1, 1])


def fdr(tool_answers, k):
    first_k = set(sorted(tool_answers[0], key=lambda _tool_answer: _tool_answer[2])[:k])
    wrong_answers = 0
    for tool_answer in first_k:
        if tool_answer.count('decoy') >= 1:
            wrong_answers += 1
    return wrong_answers / k
