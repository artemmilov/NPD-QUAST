import os
import shutil

from abstract.abstract_tool import AbstractTool

from subprocess import run, CalledProcessError


class SiriusTool(AbstractTool):
    _spectra_format = 'mgf'
    _database_format = 'txt'
    _tool_name = 'Sirius'

    def _init_tool(self, abs_folder):
        super()._init_tool(abs_folder)
        os.mkdir(
            os.path.join(abs_folder, 'temp', 'tool', 'cur_results'),
        )

    def _run_tool(self, abs_folder, specification=None):
        path_to_spectres = os.path.join(abs_folder, 'temp', 'spectres')
        path_to_results = os.path.join(abs_folder, 'temp', 'tool', 'cur_results')
        path_to_database = os.path.join(abs_folder, 'temp', 'database.txt')
        if os.path.isdir(path_to_results):
            shutil.rmtree(path_to_results)
            os.mkdir(path_to_results)

        for spectra in os.listdir(path_to_spectres):
            try:
                run(
                    '''export PATH=\"/home/artem/Programming/bioinformatics/sirius/bin/:$PATH\"; \
sirius \
-i \"{0}\" \
-o \"{1}\" \
config \
--AlgorithmProfile qtof \
--IsotopeMs2Settings IGNORE \
--MS2MassDeviation.allowedMassDeviation "10.0ppm (0.01 Da)" \
--NumberOfCandidatesPerIon 1 \
--Timeout.secondsPerTree 1000 \
--NumberOfCandidates 10 \
--FormulaSettings.enforced HCNOPS \
--Timeout.secondsPerInstance 0 \
--AdductSettings.detectable "[[M+H]+]" \
--StructureSearchDB \"{2}\" \
--AdductSettings.fallback "[[M+H]+]" \
--FormulaResultThreshold true \
--RecomputeResults true \
formula \
structure'''.format(
                        os.path.join(path_to_spectres, spectra),
                        os.path.join(path_to_results, spectra.split('.')[0]),
                        path_to_database,
                    ),
                    shell=True,
                )
            except CalledProcessError:
                pass

    def _parse_output(self, abs_folder, challenge_name):
        pass
