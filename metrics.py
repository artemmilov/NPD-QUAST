def percent_true_best(true_answers, tool_answers, best=None):
    total = len(tool_answers.keys())
    if best is None:
        best = total
    correct_matches = 0
    for scan in sorted(tool_answers, key=lambda scn: tool_answers[scn][1], reverse=True)[:best]:
        correct_matches += (tool_answers[scan][0] == true_answers[scan])
    return float(correct_matches) / min(best, total) * 100


def _abstract_medal_score(true_answers, tool_answers, medals):
    score = 0
    for scan, true_inchi in true_answers.items():
        candidate_inches = list(
            map(
                lambda scn: tool_answers[scn][0],
                sorted(
                    filter(
                        lambda tool_scan: tool_scan == scan,
                        tool_answers,
                    ),
                    key=lambda scn: tool_answers[scn][1],
                    reverse=True,
                ),
            ),
        )
        if true_inchi in candidate_inches:
            pos = candidate_inches.index(true_inchi)
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
