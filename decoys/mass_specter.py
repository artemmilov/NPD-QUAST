def _is_peak(s):
    if len(s.split('\t')) != 2:
        return False
    try:
        float(s.split('\t')[0])
    except ValueError:
        return False
    try:
        float(s.split('\t')[1])
    except ValueError:
        return False
    return True


class MassSpecter:
    peaks = []

    def __init__(self, text=None, file=None, peaks=None):
        if int(text is not None) + int(file is not None) + int(peaks is not None) > 1:
            raise AttributeError('Too many params')
        elif int(text is not None) + int(file is not None) + int(peaks is not None) == 0:
            raise AttributeError('You should use one of the params: text, file or peaks')

        self.peaks = []
        _text = ''
        if peaks is not None:
            self.peaks = peaks
            return
        elif file is not None:
            with open(file) as f:
                _text = f.read()
        elif text is not None:
            _text = text

        for line in _text.split('\n'):
            if _is_peak(line):
                self.peaks.append([float(line.split('\t')[0]), float(line.split('\t')[1])])

    def write_to_file(self, folder):
        with open(folder, 'w') as f:
            f.write('''BEGIN IONS
MSLEVEL=2
PEPMASS=150.0913
CHARGE=1+
SCANS=1\n''')
            for peak in self.peaks:
                f.write('{}\t{}\n'.format(peak[0], peak[1]))
            f.write('END IONS\n')
