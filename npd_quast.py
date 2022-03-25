import os

from tools.magma_tool import MagmaTool
from tools.NPDTools.dereplicator_tool import DereplicatorTool
from tools.NPDTools.dereplicator_plus_tool import DereplicatorPlusTool
from tools.sirius_tool import SiriusTool
from npd_quast_folder import NPDQuastFolder


def main(default_input=None):
    if default_input is None:
        args = input().split()
    else:
        args = default_input.split()
    command = args[0]
    if command == 'tool_report':
        tool_name, work_folder = args[1:]
        if tool_name == 'MAGMa+':
            tool = MagmaTool()
        elif tool_name == 'Dereplicator':
            tool = DereplicatorTool()
        elif tool_name == 'Dereplicator+':
            tool = DereplicatorPlusTool()
        elif tool_name == 'Sirius':
            tool = SiriusTool()
        else:
            print('Incorrect tool name!')
            return
        if not os.path.isdir(work_folder):
            print('There`s no directory: {}'.format(work_folder))
            return
        folder = NPDQuastFolder(work_folder)
        folder.make_tool_report(tool)
        print(
            '{0} report has been added!'.format(tool.name()),
        )
    elif command == 'total_report':
        work_folder, = args[1:]
        if not os.path.isdir(work_folder):
            print('There`s no directory: {}'.format(work_folder))
            return
        folder = NPDQuastFolder(work_folder)
        folder.make_total_report()
        print(
            'Total report has been added!',
        )


if __name__ == '__main__':
    for i in [0, 1, 2]:
        main(
            'tool_report {} sample/short'.format(
                'MAGMa+' * (i == 0) +
                'Dereplicator+' * (i == 1) +
                'Sirius' * (i == 2),
            ),
        )
    # tool_report Dereplicator+ sample/short
