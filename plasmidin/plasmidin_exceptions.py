class AmbiguousCutError(Exception):
    """
    Exception raised when a Bio.Restriction enzyme has an ambiguous cut site.
    """
    def __init__(self, enzyme, message = 'The enzyme has an ambiguous cut site'):
        self.enzyme = enzyme
        self.message = message
        super().__init__(self.message)

class CompatibleEndsError(Exception):
    """
    Exception raised when a Bio.Restriction enzyme has incompatible cut sites.
    """
    def __init__(self, message = 'The enzyme has an incompatible cut sites'):
        self.message = message
        super().__init__(self.message)