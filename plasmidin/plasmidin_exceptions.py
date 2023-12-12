class AmbiguousCutError(Exception):
    """
    Exception raised when a Bio.Restriction enzyme has an ambiguous cut site.
    """
    def __init__(self, enzyme, message = 'The enzyme has an ambiguous cut site'):
        self.enzyme = enzyme
        self.message = message
        super().__init__(self.message)