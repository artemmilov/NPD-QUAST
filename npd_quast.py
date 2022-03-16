import os

from converter_from_MAGMa_plus.magma_tool import MagmaTool
from converter_from_NPDTools.dereplicator_tool import DereplicatorTool
from converter_from_NPDTools.dereplicator_plus_tool import DereplicatorPlusTool
from converter_from_Sirius.sirius_tool import SiriusTool
from npd_quast_folder import NPDQuastFolder


def main():
    args = input().split()
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
    main()
    # tool_report Dereplicator+ sample/short
