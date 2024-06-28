"""ReadLibrary for paired-end or single-end read file set."""


class ReadLibrary():
    """Represent sequence library typically a single or paired-end FASTQ read files."""

    is_paired = False                   # Paired library status
    file_set = [None, None]              # set of files usually just 1 or 2 to store
    name = ""                           # Name of the individual strain (eg ReadGroup)

    def __init__(self, fileset: list[str], ind_name: str):
        """Initialize the ReadLibrary."""
        if fileset is not None:
            if len(fileset) > 2:
                raise ValueError("Fileset needs to be either one or two files provided")
            self.is_paired = len(fileset) == 2
            self.file_set = fileset
            self.name = ind_name

    def left(self) -> str:
        """Shortcut to get left read pair."""
        return self.file_set[0]

    def right(self) -> str:
        """Shortcut to get right read pair."""
        return self.file_set[1]
