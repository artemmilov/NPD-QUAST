def percent_true_best(true_answers, tool_answers, best=None):
    total = len(tool_answers.keys())
    if best is None:
        best = total
    correct_matches = 0
    for scan in sorted(tool_answers, key=lambda scn: tool_answers[scn][1], reverse=True)[:best]:
        correct_matches += (tool_answers[scan][0] == true_answers[scan])
    return float(correct_matches) / min(best, total) * 100


def main():
    answers_folder = input()
    true_answers = {}
    tool_answers = {}
    with open(answers_folder + '/true_answers.txt') as true_answers_data:
        for true_answer in true_answers_data.read().split('\n'):
            if true_answer != '':
                true_answers[int(true_answer.split('\t')[0])] = true_answer.split('\t')[1]
    with open(answers_folder + '/tool_answers.txt') as tool_answers_data:
        for tool_answer in tool_answers_data.read().split('\n'):
            if tool_answer != '':
                tool_answers[int(tool_answer.split('\t')[0])] = (tool_answer.split('\t')[1], tool_answer.split('\t')[2])
    for best in [10, 100, 1000, None]:
        print('{0:0.2f}% correct in best {1}.'.format(percent_true_best(true_answers, tool_answers, best=best), best))


if __name__ == '__main__':
    main()  # converter_from_NPDTools/test_main_input
