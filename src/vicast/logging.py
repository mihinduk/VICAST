"""
Logging configuration for VICAST.

This module provides standardized logging setup for VICAST tools.
It supports console output with colors and optional file logging.
"""

import logging
import sys
from pathlib import Path
from typing import Optional


# ANSI color codes for terminal output
class Colors:
    """ANSI color codes for terminal output."""
    RESET = "\033[0m"
    RED = "\033[0;31m"
    GREEN = "\033[0;32m"
    YELLOW = "\033[1;33m"
    BLUE = "\033[0;34m"
    MAGENTA = "\033[0;35m"
    CYAN = "\033[0;36m"
    GRAY = "\033[0;90m"


class ColoredFormatter(logging.Formatter):
    """
    Custom formatter with colored output for different log levels.
    """

    LEVEL_COLORS = {
        logging.DEBUG: Colors.GRAY,
        logging.INFO: Colors.GREEN,
        logging.WARNING: Colors.YELLOW,
        logging.ERROR: Colors.RED,
        logging.CRITICAL: Colors.RED,
    }

    def __init__(self, fmt: Optional[str] = None, use_colors: bool = True):
        super().__init__(fmt or "%(message)s")
        self.use_colors = use_colors and sys.stderr.isatty()

    def format(self, record: logging.LogRecord) -> str:
        # Add color to level name
        if self.use_colors:
            color = self.LEVEL_COLORS.get(record.levelno, Colors.RESET)
            record.levelname = f"{color}[{record.levelname}]{Colors.RESET}"
        else:
            record.levelname = f"[{record.levelname}]"

        return super().format(record)


def setup_logging(
    name: str = "vicast",
    level: int = logging.INFO,
    log_file: Optional[Path] = None,
    use_colors: bool = True,
    verbose: bool = False,
) -> logging.Logger:
    """
    Set up logging for VICAST tools.

    Args:
        name: Logger name (default: "vicast")
        level: Logging level (default: INFO)
        log_file: Optional path to log file
        use_colors: Use colored output for console (default: True)
        verbose: Enable verbose/debug output (default: False)

    Returns:
        Configured logger instance
    """
    if verbose:
        level = logging.DEBUG

    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Clear existing handlers
    logger.handlers.clear()

    # Console handler with colors
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(level)
    console_fmt = "%(levelname)s %(message)s"
    console_handler.setFormatter(ColoredFormatter(console_fmt, use_colors=use_colors))
    logger.addHandler(console_handler)

    # File handler (plain text, no colors)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_fmt = "%(asctime)s %(levelname)s [%(name)s] %(message)s"
        file_handler.setFormatter(logging.Formatter(file_fmt))
        logger.addHandler(file_handler)

    return logger


def get_logger(name: str = "vicast") -> logging.Logger:
    """
    Get or create a logger for VICAST tools.

    Args:
        name: Logger name (default: "vicast")

    Returns:
        Logger instance
    """
    logger = logging.getLogger(name)

    # If no handlers, set up with defaults
    if not logger.handlers:
        setup_logging(name)

    return logger


# Convenience function for scripts
def log_info(message: str, logger_name: str = "vicast") -> None:
    """Log an info message."""
    get_logger(logger_name).info(message)


def log_warning(message: str, logger_name: str = "vicast") -> None:
    """Log a warning message."""
    get_logger(logger_name).warning(message)


def log_error(message: str, logger_name: str = "vicast") -> None:
    """Log an error message."""
    get_logger(logger_name).error(message)


def log_debug(message: str, logger_name: str = "vicast") -> None:
    """Log a debug message."""
    get_logger(logger_name).debug(message)


# Create default logger on import
_default_logger = setup_logging()
