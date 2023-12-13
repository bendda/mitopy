import logging
import sys
from shutil import which
import subprocess


class Executable:
    def __init__(self, exec_path, verbosity=False):
        self.exec_path = exec_path
        self.verbosity = verbosity

    def _check_executable(self):
        if not which(self.exec_path):
            logging.error(f"Executable {self.exec_path} not found.")
            return False
        return True

    def _build_cmd(self, *args, subcommand=None, output_redirect=None, **kwargs) -> str:
        """Build a command with provided positional arguments and other options.

        Returns:
            str: Created command.
        """
        cmd = f"{self.exec_path} {subcommand} " if subcommand else f"{self.exec_path} "

        for arg, value in kwargs.items():
            if isinstance(value, bool) and value:
                cmd += f"{arg} "
            elif isinstance(value, list):
                cmd += " ".join([f"{arg} {opt}" for opt in value]) + " "
            else:
                cmd += f"{arg} {value} "

        cmd += " ".join(args)
        if output_redirect:
            cmd += f" > {output_redirect}"
        return cmd

    def run(
        self,
        *args,
        subcommand=None,
        redirect_out=None,
        verbose=False,
        **kwargs,
    ) -> bool:
        if self._check_executable():
            # Build command
            cmd = self._build_cmd(
                *args,
                subcommand=subcommand,
                output_redirect=redirect_out,
                **kwargs,
            )

            # Execute command
            try:
                logging.info(f"Running command: {cmd}")
                subprocess.run(
                    cmd,
                    check=True,
                    capture_output=(not self.verbosity),
                    shell=True,
                    text=True,
                )
            except subprocess.CalledProcessError as e:
                logging.error(f"{cmd} failed with following error: \n\n {e.stderr}")
                sys.exit(1)
        else:
            sys.exit(1)
