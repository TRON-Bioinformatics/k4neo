from logzero import logger
import os
from snakemake import snakemake


class Configuration:
    """
    neoKant configuration object to manage annotation pipeline and others
    """

    def __init__(self, config_file, verbose=False):
        self.config = snakemake.load_configfile(config_file)
        if verbose:
            self.log_configuration()

    def log_configuration(self):
        logger.info("Configuration")
        for k, v in self.config.items():
            logger.info("{}={}".format(k, v))
