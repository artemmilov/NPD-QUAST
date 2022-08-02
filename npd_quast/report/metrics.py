import numpy as np
from numpy import mean, median, quantile
from rdkit import Chem, DataStructs

ROUND = 2


def top_x(true_answers, tool_answers, x=None):
    all_answers = []
    for challenge, challenge_answers in tool_answers.items():
        for challenge_answer in challenge_answers:
            all_answers.append(
                (
                    challenge,
                    challenge_answer[0],
                    challenge_answer[2],
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
        lambda ans: ans[1] == true_answers[ans[0]][0],
        all_answers,
    )))


def _sort_by_score(tool_answers, scan):
    sorted_answers = sorted(
        [
            (tool_answer[0], tool_answer[1], tool_answer[2],)
            for tool_answer in tool_answers[scan]
        ],
        key=lambda ans: ans[2],
    )
    result = []
    last_scan = None
    for sorted_answer in sorted_answers:
        if sorted_answer[2] == last_scan:
            result[-1].append((sorted_answer[0], sorted_answer[1]))
        else:
            result.append([(sorted_answer[0], sorted_answer[1])])
            last_scan = sorted_answer[2]
    return result


def _get_ranks(true_answers, tool_answers, default_rank=None):
    ranks = []
    for scan, [true_inchi, _] in true_answers.items():
        sorted_answers = _sort_by_score(tool_answers, scan)
        found = False
        for i, set_answers in enumerate(sorted_answers):
            if true_inchi in map(lambda ans: ans[0], set_answers):
                before = sum(map(len, sorted_answers[:i]))
                ans_rank = mean(
                    list(
                        map(
                            lambda x: before + x,
                            range(0, len(set_answers)),
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
    return round(
        float(mean(_get_ranks(true_answers, tool_answers, default_rank))),
        ROUND,
    )


def median_rank(true_answers, tool_answers, default_rank=None):
    return round(
        float(
            median(_get_ranks(true_answers, tool_answers, default_rank))
        ),
        ROUND,
    )


def _get_rrps(true_answers, tool_answers):
    rrps = []
    for scan, [true_inchi, _] in true_answers.items():
        sorted_answers = _sort_by_score(tool_answers, scan)
        for i, set_answers in enumerate(sorted_answers):
            if true_inchi in map(lambda ans: ans[0], set_answers):
                before = sum(map(len, sorted_answers[:i]))
                after = sum(map(len, sorted_answers[i + 1:]))
                total = sum(list(map(len, sorted_answers)))
                if total > 1:
                    rrps.append(1 / 2 * (1 - (before - after) / (total - 1)))
                else:
                    rrps.append(1 / 2)
    return rrps


def mean_rrp(true_answers, tool_answers):
    return round(
        float(mean(_get_rrps(true_answers, tool_answers))),
        ROUND,
    )


def median_rrp(true_answers, tool_answers):
    return round(
        float(median(_get_rrps(true_answers, tool_answers))),
        ROUND,
    )


def _get_weighted_rrps(true_answers, tool_answers):
    weighted_rrps = []
    for scan, [true_inchi, _] in true_answers.items():
        sorted_answers = sorted(
            [
                [tool_answer[0], tool_answer[2]]
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
    return round(
        float(
            mean(_get_weighted_rrps(true_answers, tool_answers))
        ),
        ROUND,
    )


def median_weighted_rrp(true_answers, tool_answers):
    return round(
        float(
            median(_get_weighted_rrps(true_answers, tool_answers))
        ),
        ROUND,
    )


def k_quantile(true_answers, tool_answers, k=50):
    positions = []
    for scan, [true_inchi, _] in true_answers.items():
        sorted_answers = _sort_by_score(tool_answers, scan)
        for i, set_answers in enumerate(sorted_answers):
            if true_inchi in map(lambda ans: ans[0], set_answers):
                before = sum(map(len, sorted_answers[:i]))
                ans_rank = mean(
                    list(
                        map(
                            lambda x: before + x,
                            range(0, len(set_answers)),
                        ),
                    ),
                )
                positions.append(ans_rank)
    if len(positions) != 0:
        return round(quantile(positions, k / 100), ROUND)
    return 0


def _abstract_medal_score(true_answers, tool_answers_dict, medals):
    scores_dict = {tool_title: 0 for tool_title in tool_answers_dict}
    for scan, [true_inchi, _] in true_answers.items():
        positions_dict = {}
        for tool_title, tool_answers in tool_answers_dict.items():
            sorted_tool_answers = sorted(
                tool_answers[scan],
                key=lambda ans: ans[2],
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


def _list_similarity_top_x(true_answers, tool_answers, x):
    similarity = {}
    for scan, [true_inchi, true_smiles] in true_answers.items():
        m = Chem.MolFromSmiles(true_smiles)
        true_fingerprint = Chem.RDKFingerprint(m)
        sorted_answers = _sort_by_score(tool_answers, scan)
        cur_similarity = 0
        for i, set_answers in enumerate(sorted_answers):
            before = sum(map(len, sorted_answers[:i]))
            ans_rank = mean(
                list(
                    map(
                        lambda x: before + x,
                        range(0, len(set_answers)),
                    ),
                ),
            )
            if (x is not None) and (ans_rank > x):
                break
            for answer in set_answers:
                m = Chem.MolFromSmiles(answer[1])
                cur_fingerprint = Chem.RDKFingerprint(m)
                cur_similarity = max(
                    cur_similarity,
                    DataStructs.FingerprintSimilarity(
                        true_fingerprint,
                        cur_fingerprint,
                    ),
                )
            similarity[scan] = cur_similarity
    return similarity


def mean_similarity_top_x(true_answers, tool_answers, x=None):
    return np.mean(list(_list_similarity_top_x(true_answers, tool_answers, x).values()))


def median_similarity_top_x(true_answers, tool_answers, x=None):
    return np.median(list(_list_similarity_top_x(true_answers, tool_answers, x).values()))
