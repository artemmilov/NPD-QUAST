import math
import os

from decoys.mass_specter import MassSpecter
import random
from collections import defaultdict


def handle_naive_method(mass_spectra, decoys_count=1, number_peaks_like_target=False, peaks_from_same_pepmass=False):
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
    result_spectra = []
    for _ in range(decoys_count):
        result_peaks = []
        for _ in range(k):
            result_peaks.append(random.choice(all_spectres))
        result_spectra.append(MassSpecter(peaks=result_peaks))
    return result_spectra


def handle_naive_method_extra(mass_spectra, decoys_count=1, number_peaks_like_target=False,
                              peaks_count=100):
    peak_freq = defaultdict(int)
    for mass_specter in mass_spectra:
        for peak in mass_specter.peaks:
            peak_freq[str(int(100 * peak[0]))] += 1

    lower_bound = min(sorted(peak_freq.values(), reverse=True)[:peaks_count])
    for peak, freq in peak_freq.items():
        if freq < lower_bound:
            peak_freq[peak] = 0

    total = sum(peak_freq.values())
    for peak in peak_freq:
        peak_freq[peak] /= total
    result_spectra = []

    for i in range(decoys_count):
        if int(i / decoys_count * 49) > int((i - 1) / decoys_count * 49):
            print('#' * (int(i / decoys_count * 49) - int((i - 1) / decoys_count * 49)), end='', flush=True)
        target_mass_specter = random.choice(mass_spectra)
        target_length = len(target_mass_specter.peaks)
        k = target_length
        # if not number_peaks_like_target:
        #     k = random.choice(list(map(lambda ms: len(ms.peaks), mass_spectra)))
        # else:
        #     k = target_length
        result_peaks = []
        for _ in range(k):
            peak = random.choices(list(peak_freq.keys()), weights=list(peak_freq.values()))
            result_peaks.append([peak, 1.0])
        result_spectra.append(MassSpecter(peaks=result_peaks))
    return result_spectra


def handle_spectrum_based_method(mass_spectra, decoys_count=1,
                                 number_peaks_like_target=False, filter_by_last_peak=False):
    result_spectra = []
    for _ in range(decoys_count):
        target_mass_specter = random.choice(mass_spectra)
        target_pepmass = target_mass_specter.pepmass
        target_length = len(target_mass_specter.peaks)

        if not number_peaks_like_target:
            k = random.choice(list(map(lambda ms: len(ms.peaks), mass_spectra)))
        else:
            k = target_length

        result_peaks = [random.choice(target_mass_specter.peaks)]
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
        result_spectra.append(MassSpecter(peaks=result_peaks))
    return result_spectra


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
                                                                         'decoy_{}.mgf'.format(i)))


def main():
    npd_quast_folder = '../../NPD-QUAST-test/snakemake_test/test_decoys'
    make_decoys(npd_quast_folder)


if __name__ == '__main__':
    main()
