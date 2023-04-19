import os

from decoys.mass_specter import MassSpecter
import random


def handle_naive_method(mass_spectra, number_peaks_like_target=False, peaks_from_same_pepmass=False):
    target_mass_specter = random.choice(mass_spectra)
    target_pepmass = target_mass_specter.pepmass
    target_length = len(target_mass_specter.peaks)

    all_spectres = []
    for mass_specter in mass_spectra:
        if (not peaks_from_same_pepmass) or (mass_specter.pepmass == target_pepmass):
            for peak in mass_specter.peaks:
                all_spectres.append(peak)
    if not number_peaks_like_target:
        k = random.choice(list(map(lambda ms: len(ms.peaks), mass_spectra)))
    else:
        k = target_length
    result_specter = []
    for _ in range(k):
        result_specter.append(random.choice(all_spectres))
    return MassSpecter(peaks=result_specter)


def handle_spectrum_based_method(mass_spectra, number_peaks_like_target=False, filter_by_last_peak=False):
    target_mass_specter = random.choice(mass_spectra)
    target_pepmass = target_mass_specter.pepmass
    target_length = len(target_mass_specter.peaks)
    result_peaks = [random.choice(target_mass_specter.peaks)]

    if not number_peaks_like_target:
        k = random.choice(list(map(lambda ms: len(ms.peaks), mass_spectra)))
    else:
        k = target_length

    for _ in range(k - 1):
        if filter_by_last_peak:
            search_peak_mass = result_peaks[-1][0]
        else:
            search_peak_mass = random.choice(result_peaks)[0]
        suitable_peaks = []
        for mass_specter in mass_spectra:
            if (search_peak_mass in set(map(lambda pk: pk[0], mass_specter.peaks))) and \
                    (mass_specter.pepmass == target_pepmass):
                suitable_peaks += mass_specter.peaks
        result_peaks.append(random.choice(suitable_peaks))
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
