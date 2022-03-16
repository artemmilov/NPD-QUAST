from tools.magma_tool import MagmaTool
from tools.NPDTools.dereplicator_tool import DereplicatorTool
from tools.NPDTools.dereplicator_plus_tool import DereplicatorPlusTool
from tools.sirius_tool import SiriusTool
from npd_quast_folder import NPDQuastFolder


def main():
    folder = None
    print("""Which folder do you want work with?
Short:\t's'
Medium:\t'm'
Large:\t'l'
Exit:\t'e'""")
    ok = False
    while not ok:
        a = input()
        ok = True
        if a == 's':
            folder = NPDQuastFolder('short')
        elif a == 'm':
            folder = NPDQuastFolder('medium')
        elif a == 'l':
            folder = NPDQuastFolder('large')
        elif a == 'e':
            return
        else:
            ok = False
            print('Sorry?!')
    first = True
    print("""Which tool do you want report?
MAGMa+:\t\t\t0
Dereplicator:\t1
Dereplicator+:\t2
Sirius:\t\t\t3
exit:\t\t\t'e'""")
    while True:
        a = input()
        if a == '0':
            tool = MagmaTool()
        elif a == '1':
            tool = DereplicatorTool()
        elif a == '2':
            tool = DereplicatorPlusTool()
        elif a == '3':
            tool = SiriusTool()
        elif (a == 'g') and (not first):
            folder.make_total_report()
            print("""Global report has been added!\n
Which tool do you want report further?
MAGMa+:\t0
Dereplicator:\t1
Dereplicator+:\t2
Sirius:\t3
Global:\t'g'
Exit:\t'e'""")
            continue
        elif a == 'e':
            return
        else:
            print('Sorry?!')
            continue
        folder.make_tool_report(tool)
        print("""{0} report has been added!\n
Which tool do you want report further?
MAGMa+:\t0
Dereplicator:\t1
Dereplicator+:\t2
Sirius:\t3
Global:\t'g'
Exit:\t'e'""".format(tool.name()))
        first = False


if __name__ == '__main__':
    main()
