from rdkit import Chem


class DereplicatorSpectreIdentification:
    def __init__(self, spec_file, local_spec_idx, scan, local_peptide_idx, name, score, p_value):
        self.spec_file = spec_file
        self.local_spec_idx = local_spec_idx
        self.scan = scan
        self.local_peptide_idx = local_peptide_idx
        self.name = name
        self.score = score
        self.p_value = p_value


class TrueSpectreIdentification:
    def __init__(self, local_spec_idx, scan, local_peptide_idx, name, score):
        self.local_spec_idx = local_spec_idx
        self.scan = scan
        self.local_peptide_idx = local_peptide_idx
        self.name = name
        self.score = score


def initialize_dereplicator_specs(all_matches_file):
    spectre_identifications = []
    with open(all_matches_file) as all_matches:
        for match in map(lambda s: s.split('\t'), all_matches.readlines()[1:]):
            spectre_identifications.append(DereplicatorSpectreIdentification(
                match[0], match[1], match[2], match[3], match[4], match[5], match[6]
            ))
    return spectre_identifications


def initialize_true_specs(library_file, smiles_info_file):
    spectre_identifications = []
    with open(library_file) as library, \
            open(smiles_info_file) as smiles_info:
        all_smiles = smiles_info.read().splitlines()
        for spectra in library.read().split('\n\n'):
            local_spec_idx, scan, local_peptide_idx = 0, 0, 0
            score = 0
            smiles, name, mass = None, None, None
            scan = 0
            for spectra_line in spectra.splitlines():
                if '=' in spectra_line:
                    spectra_line_split = spectra_line.split('=', maxsplit=1)
                    left, right = spectra_line_split[0], spectra_line_split[1]
                    if (left == 'SMILES') and ('.' not in right):
                        smiles = right
                    elif (left == 'NAME') and (name != ''):
                        name = right.replace(' ', '_')
                    elif left == 'PEPMASS':
                        mass = float(right)
            if (name is not None) and \
                    (mass is not None):
                ok = True
                try:
                    m = Chem.MolFromSmiles(smiles)
                    m = Chem.AddHs(m)
                except Exception:
                    ok = False
                if ok:
                    local_peptide_idx = all_smiles.index(smiles)
                    spectre_identifications.append(TrueSpectreIdentification(
                        local_spec_idx, scan, local_peptide_idx, name, score
                    ))
                    local_spec_idx += 1
        scan += 1
    return spectre_identifications


def _is_identification_correct(dereplicator_identification, true_identifications):
    for true_identification in true_identifications:
        if dereplicator_identification.local_peptide_idx == true_identification.local_peptide_idx:
            return True
    return False


def percent_true_top(spectre_identifications, true_identifications, best=None):
    if best is None:
        best = len(spectre_identifications)
    true_answers = 0
    for spectre_identification in sorted(spectre_identifications, key=lambda si: si.p_value, reverse=True)[:best]:
        true_answers += _is_identification_correct(spectre_identification, true_identifications)
    return float(true_answers) / best * 10


def main():
    try:
        library_address, all_matches_address, smiles_address = input().split()
    except (ValueError, EOFError):
        print('Incorrect input!')
        return
    true_identifications = initialize_true_specs(library_address, smiles_address)
    dereplicator_identifications = initialize_dereplicator_specs(all_matches_address)
    print('{} {} {} {}'.format(
        percent_true_top(dereplicator_identifications, true_identifications),
        percent_true_top(dereplicator_identifications, true_identifications, best=10),
        percent_true_top(dereplicator_identifications, true_identifications, best=100),
        percent_true_top(dereplicator_identifications, true_identifications, best=1000)
    ))


if __name__ == '__main__':
    main()  # converter_from_NPDTools/input_sample/input_library.mgf converter_from_NPDTools/test/all_matches.tsv converter_from_NPDTools/output_sample/database/smiles.info
