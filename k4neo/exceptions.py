class K4neoException(Exception):
    pass


class K4neoInputException(K4neoException):
    pass


class K4neoPipelineException(K4neoException):
    pass


class K4neoParsingException(K4neoException):
    pass


class K4neoParsingStudyException(K4neoParsingException):
    pass
