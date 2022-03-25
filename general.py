from collections import defaultdict


def parse_from_mgf(s):
    res = []
    cur = ''
    in_brackets = False
    for a in s:
        if (a == ',') and (not in_brackets):
            res.append(cur)
            cur = ''
        elif a == '\"':
            in_brackets = not in_brackets
        else:
            cur += a
    if cur != '':
        res.append(cur)
    return res


def get_name(conn, user_id):
    c = conn.cursor()
    c.execute("SELECT name FROM userList WHERE userID = ?", (user_id,))
    result = c.fetchone()
    if result:
        return result[0]


def parse_true_answers(true_answers_data_file):
    true_answers = {}
    with open(true_answers_data_file) as true_answers_data:
        for true_answer in true_answers_data.read().split('\n'):
            if true_answer != '':
                true_answers[true_answer.split('\t')[0]] = \
                    true_answer.split('\t')[1]
    return true_answers


def parse_tool_answers(tool_answers_data_file):
    tool_answers = defaultdict(list)
    with open(tool_answers_data_file) as tool_answers_data:
        for tool_answer in tool_answers_data.read().split('\n'):
            if tool_answer != '':
                tool_answers[tool_answer.split('\t')[0]].append(
                    (
                        tool_answer.split('\t')[1],
                        tool_answer.split('\t')[2],
                    ),
                )
    return tool_answers
