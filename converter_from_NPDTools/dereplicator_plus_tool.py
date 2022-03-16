import os
from subprocess import run

from converter_from_NPDTools.abstract_npd_tool import AbstractNpdTool


class DereplicatorPlusTool(AbstractNpdTool):
    _tool_name = 'Dereplicator_plus'

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
                '--pass-to-dereplicate',
                '--num_hits_to_report 100',
                '--min-score',
                '1',
            ],
        )
