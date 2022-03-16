from general import parse_from_mgf


def create_true_answers(web_results_file, true_answers_file):
    with open(web_results_file) as web_results, \
         open(true_answers_file, 'w') as true_answers:
        for row in web_results.readlines()[1:]:
            challenge = result = parse_from_mgf(row)[1]
            answer = parse_from_mgf(row)[9].split('-')[0]
            true_answers.write(
                '{0}${1}\t{2}\n'.format(
                    challenge,
                    result,
                    answer
                ),
            )


def main():
    web_results_file, true_answers_file = input().split()
    create_true_answers(web_results_file, true_answers_file)


if __name__ == '__main__':
    main()
