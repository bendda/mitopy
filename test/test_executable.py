from mitopy.executable import Executable
import logging
import pytest


def test_check_executable(caplog):
    """Unit test for check_dependencies method."""
    caplog.set_level(logging.INFO)

    exec = Executable(exec_path="ls")
    status = exec._check_executable()
    assert status == True

    exec = Executable(exec_path="non_existent_exec")
    status = exec._check_executable()
    assert f"Executable non_existent_exec not found." in caplog.text


def test_build_cmd(caplog):
    """Unit test for build_cmd method."""
    exec = Executable(exec_path="ls")

    assert (
        exec._build_cmd(subcommand="mem", output_redirect="output.txt", **{"-l": True})
        == "ls mem -l  > output.txt"
    )
    assert (
        exec._build_cmd("file", subcommand="mem", output_redirect="output.txt")
        == "ls mem file > output.txt"
    )
    assert (
        exec._build_cmd(
            "file", subcommand="mem", output_redirect="output.txt", **{"-l": True}
        )
        == "ls mem -l file > output.txt"
    )
    assert exec._build_cmd("file", subcommand=None, output_redirect=None) == "ls file"
    assert (
        exec._build_cmd(
            "file",
            subcommand=None,
            output_redirect=None,
            **{"--test": "test"},
        )
        == "ls --test test file"
    )


def test_run(caplog):
    """Unit test for Executable run method."""
    caplog.set_level(logging.INFO)

    exec = Executable(exec_path="ls")
    assert exec.exec_path == "ls"
    params = {"-l": True, "-a": True}
    exec.run(**params)
    assert "Running command: ls -l -a" in caplog.text
