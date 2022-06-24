import os
import json

import npd_quast.tools as tools
import npd_quast.report
import npd_quast.npd_quast_folder as npd_quast_folder
import npd_quast.general as general


def run_n_report(options):
    if options.tool in tools.SUPPORTED_TOOLS.keys():
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
                folder = npd_quast_folder.NPDQuastFolder(options.folder)
            except AttributeError as e:
                print(e)
                return
            except NotADirectoryError as e:
                print(e)
                return
            folder.make_tool_report(tool, options.report_name, specification)
            print(
                '{0} report has been reported in {1}!'.format(
                    options.tool,
                    options.report_name,
                ),
            )


def compile_reports(options):
    try:
        npd_quast_folder.NPDQuastFolder(options.folder)
    except AttributeError as e:
        print(e)
        return
    except NotADirectoryError as e:
        print(e)
        return
    if os.path.isdir(options.folder):
        abs_folder = os.path.abspath(
            os.path.join(
                options.folder,
            )
        )
        true_answers = general.parse_true_answers(
            os.path.join(
                abs_folder,
                'true_answers.txt',
            )
        )
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
            true_answers,
            tool_answers_dict,
        )
        print('All reports has been compiled!')
