class NeoKantException(Exception):
    pass


class NeoKantPipelineException(NeoKantException):
    pass


class NeoKantParsingException(NeoKantException):
    pass


class NeoKantParsingStudyException(NeoKantParsingException):
    pass
