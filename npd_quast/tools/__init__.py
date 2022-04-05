from .magma_tool import MagmaTool
from .npdtools import DereplicatorPlusTool, DereplicatorTool
from .sirius_tool import SiriusTool

SUPPORTED_TOOLS = {
    'MAGMa+': MagmaTool,
    'Dereplicator': DereplicatorTool,
    'Dereplicator+': DereplicatorPlusTool,
    'Sirius': SiriusTool,
}
