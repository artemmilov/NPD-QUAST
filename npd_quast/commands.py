import os
import json

import npd_quast.tools as tools
import npd_quast.report
import npd_quast.npd_quast_folder as npd_quast_folder
import npd_quast.general as general


def run_n_report(options, logger):
    if options.tool in tools.SUPPORTED_TOOLS.keys():
        logger.info('Running \"{}\" tool'.format(options.tool))
        tool = tools.SUPPORTED_TOOLS[options.tool]()
        if os.path.isdir(options.folder):
            if options.config is not None:
                f = open(options.config)
            else:
                f = open(
                    os.path.abspath(
                        os.path.join(
                            os.curdir,
                            'default_configurations',
                            tool.name() + '.json',
                        )
                    )
                )
            specification = json.load(f)
            try:
                folder = npd_quast_folder.NPDQuastFolder(options.folder, logger)
            except AttributeError or NotADirectoryError:
                logger.error(' ', is_exception=True)
                return
            logger.info('Input data is ok. Started making report...')
            folder.make_tool_report(
                tool,
                options.report_name,
                specification,
                logger,
                debug=options.debug,
            )


def compile_reports(options, logger):
    try:
        npd_quast_folder.NPDQuastFolder(options.folder, logger)
    except AttributeError or NotADirectoryError:
        logger.error(' ', is_exception=True)
        return
    if os.path.isdir(options.folder):
        abs_folder = os.path.abspath(
            os.path.join(
                options.folder,
            )
        )

        true_answers_on = 'true_answers.txt' in os.listdir(options.folder)
        true_answers = None
        if true_answers_on:
            true_answers = general.parse_true_answers(
                os.path.join(
                    abs_folder,
                    'true_answers.txt',
                )
            )
        decoys_on = False
        for challenge in os.listdir(os.path.join(options.folder, 'challenges')):
            decoys_on = decoys_on or 'decoys' in os.listdir(os.path.join(options.folder, 'challenges', challenge))
        tool_answers_dict = {
            report: general.parse_tool_answers(
                os.path.join(
                    abs_folder,
                    'reports',
                    report,
                    'tool_answers.txt',
                ),
            )
            for report in os.listdir(
                os.path.join(
                    abs_folder,
                    'reports',
                )
            )
            if os.path.isdir(
                os.path.join(
                    abs_folder,
                    'reports',
                    report,
                ),
            )
        }
        npd_quast.report.write_report(
            abs_folder,
            tool_answers_dict,
            true_answers=true_answers,
            true_answers_on=true_answers_on,
            decoys_on=decoys_on
        )
        print('All reports has been compiled!')
