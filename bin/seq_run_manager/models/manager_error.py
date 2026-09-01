class ManagerError(Exception):
    """Exception raised within the manager."""

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
