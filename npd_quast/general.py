import os.path
from collections import defaultdict


class NPDQuastError(Exception):
    pass


def get_true_spectra(folder):
    true_spectra = []
    for challenge in os.listdir(os.path.join(folder, 'challenges')):
        if challenge != 'fake':
            for specter in os.listdir(os.path.join(folder, 'challenges', challenge, 'spectra')):
                true_spectra.append('{}\t{}'.format(challenge, specter))


def parse_true_answers(true_answers_data_file):
    true_answers = {}
    with open(true_answers_data_file) as true_answers_data:
        for true_answer in true_answers_data.read().split('\n'):
            if true_answer != '':
                true_answers[true_answer.split('\t')[0] + '\t' + true_answer.split('\t')[1]] = \
                    true_answer.split('\t')[2]
    return true_answers


def parse_tool_answers(tool_answers_data_file):
    tool_answers = defaultdict(list)
    with open(tool_answers_data_file) as tool_answers_data:
        for tool_answer in tool_answers_data.read().split('\n'):
            if (tool_answer != '') and (tool_answer.split('\t')[1] == 'spectra'):
                tool_answers[tool_answer.split('\t')[0] + '\t' + \
                             tool_answer.split('\t')[2]].append(
                    (
                        tool_answer.split('\t')[3],
                        float(tool_answer.split('\t')[4], )
                    ),
                )
    return tool_answers


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
