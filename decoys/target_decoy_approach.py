import os

from decoys.mass_specter import MassSpecter
import random


def handle_naive_method(mass_spectra):
    all_spectres = []
    for mass_specter in mass_spectra:
        for peak in mass_specter.peaks:
            all_spectres.append(peak)
    k = random.choice(list(map(lambda ms: len(ms.peaks), mass_spectra)))
    result_specter = []
    for _ in range(k):
        result_specter.append(random.choice(all_spectres))
    return MassSpecter(peaks=result_specter)


def handle_spectrum_based_method(mass_spectra):
    k = random.choice(list(map(lambda ms: len(ms.peaks), mass_spectra)))
    i = random.choice(range(len(mass_spectra)))
    result_peaks = [max(mass_spectra[i].peaks, key=lambda pk: pk[0])]
    result_masses = set(map(lambda pk: pk[0], result_peaks))

    for _ in range(k):
        suitable_peaks = []
        for mass_specter in mass_spectra:
            if result_masses.issubset(set(mass_specter.peaks)):
                suitable_peaks += mass_specter.peaks
        result_peaks.append(random.choice(suitable_peaks))
        result_masses = set(map(lambda pk: pk[0], result_peaks))
    return MassSpecter(peaks=result_peaks)


def make_decoys(npd_quast_folder, percent=100):
    mass_spectra = []
    for challenge in os.listdir(os.path.join(npd_quast_folder, 'challenges')):
        for specter in os.listdir(os.path.join(npd_quast_folder, 'challenges', challenge, 'spectra')):
            mass_spectra.append(MassSpecter(file=os.path.join(
                npd_quast_folder, 'challenges', challenge, 'spectra', specter)))
            print(len(MassSpecter(file=os.path.join(
                npd_quast_folder, 'challenges', challenge, 'spectra', specter)).peaks))
    for challenge in os.listdir(os.path.join(npd_quast_folder, 'challenges')):
        k = int((percent / 100) * \
                len(os.listdir(os.path.join(npd_quast_folder, 'challenges', challenge, 'spectra'))))
        for i in range(k):
            handle_naive_method(mass_spectra).write_to_file(os.path.join(npd_quast_folder,
                                                                         'challenges', challenge, 'spectra',
                                                                         'decoy_{}'.format(i)))


def main():
    npd_quast_folder = '../../NPD-QUAST-test/snakemake_test/test_decoys'
    make_decoys(npd_quast_folder)


if __name__ == '__main__':
    main()
