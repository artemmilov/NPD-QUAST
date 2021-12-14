import os


class AbstractDatabase:
    _l = None

    def __init__(self, folder):
        if len(os.listdir(folder)) == 0:
            self._empty_init(folder)
        elif len(os.listdir(folder)) != 0:
            self._not_empty_init(folder)

    def _empty_init(self, folder):
        pass

    def _not_empty_init(self, folder):
        pass

    def clean(self):
        pass

    def __len__(self):
        return self._l

    def __iter__(self):
        self._counter = 0
        return self

    def __next__(self):
        if self._counter < len(self):
            self._counter += 1
            return self[self._counter - 1]
        else:
            raise StopIteration

    def __getitem__(self, i):
        pass

    def add_mol(self, quast_mol):
        pass

    def add_database(self, other_database):
        for abstract_mol in other_database:
            self.add_mol(abstract_mol)
