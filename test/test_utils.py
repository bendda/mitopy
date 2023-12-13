from mitopy.utils import check_files_exist
import logging


def test_check_files_exist(caplog):
    """Unit test for check_files_exist function."""
    caplog.set_level(logging.INFO)

    assert check_files_exist("test/conftest.py", verbose=True) == True
    assert check_files_exist("non_existent_file.py", verbose=True) == False
    assert "File non_existent_file.py not found." in caplog.text
