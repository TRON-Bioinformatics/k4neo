import pathlib
from loguru import logger
from tqdm import tqdm


def setup_logging(log_file: pathlib.Path, verbose: bool = False):
    """
    Setup k4neo logger
    """
    logger.remove()

    # Console: Set level for CLI based on verbose flag
    console_level = "DEBUG" if verbose else "INFO"

    # Add a tqdm write as new sink. Ensures that progress bar and logs fit together
    logger.add(lambda msg: tqdm.write(msg, end=""), level=console_level, colorize=True)

    # Add sink for log file
    logger.add(
        log_file,
        level="DEBUG",
        format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
    )
    return logger