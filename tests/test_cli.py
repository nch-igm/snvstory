import os
import re
from unittest.mock import patch

import igm_churchill_cnvloh.cli as cli


@patch('os.makedirs')
def test_setup_workspace(mockdirs):
    os.environ['TMP_DIR'] = '/data'
    wrk_dir = cli.setup_workspace()
    mockdirs.assert_called_once()
    assert re.fullmatch(os.environ['TMP_DIR'] + r'/[^/]+', wrk_dir)
