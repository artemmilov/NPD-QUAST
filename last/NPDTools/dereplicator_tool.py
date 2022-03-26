from subprocess import run

from last.NPDTools.abstract_npd_tool import AbstractNpdTool


class DereplicatorTool(AbstractNpdTool):
    _tool_name = 'Dereplicator'

    def _run_tool(self, abs_folder, specification=None):
        path_to_spectres, path_to_database, path_to_result =\
            super()._run_abstract_tool(abs_folder, specification)
        run(
            [
                self._location,
                path_to_spectres,
                '--db-path',
                path_to_database,
                '-o',
                path_to_result,
            ],
        )
