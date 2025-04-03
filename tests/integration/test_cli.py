"""
Integration tests for the command-line interface.
"""
import os
import tempfile
import subprocess


def run_command(command):
    """Run a command and return the output."""
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        universal_newlines=True,
    )
    stdout, stderr = process.communicate()
    return process.returncode, stdout, stderr


class TestCLI:
    """Tests for the command-line interface."""

    def test_help_command(self):
        """Test the help command."""
        returncode, stdout, stderr = run_command("python -m gfftk --help")
        assert returncode == 0
        assert "usage:" in stdout
        # The module name might be __main__.py instead of gfftk
        assert "GFFtk" in stdout

    def test_version_command(self):
        """Test the version command."""
        returncode, stdout, stderr = run_command("python -m gfftk --version")
        assert returncode == 0
        assert "gfftk" in stdout or "GFFtk" in stdout
