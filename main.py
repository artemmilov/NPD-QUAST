import os
from collections import defaultdict

from rdkit import RDLogger

from report import write_report


def _check_folders(answers_folder, report_folder):
    if len(os.listdir(answers_folder)) != 2:
        return False
    if 'true_answers.txt' not in os.listdir(answers_folder):
        return False
    if 'tool_answers.txt' not in os.listdir(answers_folder):
        return False
    if len(os.listdir(report_folder)) != 0:
        return False
    return True


def handle_answers(answers_folder):
    true_answers = {}
    tool_answers = defaultdict(list)
    with open(os.path.join(answers_folder, 'true_answers.txt')) as true_answers_data:
        for true_answer in true_answers_data.read().split('\n'):
            if true_answer != '':
                true_answers[int(true_answer.split('\t')[0])] = true_answer.split('\t')[1]
    with open(os.path.join(answers_folder, 'tool_answers.txt')) as tool_answers_data:
        for tool_answer in tool_answers_data.read().split('\n'):
            if tool_answer != '':
                tool_answers[int(tool_answer.split('\t')[0])].append(
                    (
                        tool_answer.split('\t')[1],
                        tool_answer.split('\t')[2],
                    ),
                )
    return true_answers, tool_answers


def main():
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    while True:
        try:
            answers_folder, report_folder = input().split(' ')
            if _check_folders(answers_folder, report_folder):
                break
            print('Incorrect input!')
        except (ValueError, EOFError, SyntaxError):
            print('Incorrect input!')
    true_answers, tool_answers = handle_answers(answers_folder)
    write_report(report_folder, true_answers, tool_answers)


if __name__ == '__main__':
    main()  # converter_from_NPDTools/test_main_input test_main_output      !!!
