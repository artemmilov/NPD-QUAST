def create_mass_tree(mgf_data_file, mass_tree_file):
    with open(mgf_data_file) as mgf_data:
        mzs = []
        norm_intensities = []
        mass = None
        for line in mgf_data.readlines():
            if '=' in line:
                if line.split('=')[0] == 'PEPMASS':
                    mass = float(line.split('=')[1])
            elif len(line.split('\t')) == 2:
                mzs.append(float(line.split('\t')[0]))
                norm_intensities.append(float(line.split('\t')[1]))
        normalizer = max(norm_intensities)
        record = '{0}: 100 ('.format(mass)
        record += ', '.join(
            '{0}: {1}'.format(
                mzs[i],
                norm_intensities[i] * 100 / normalizer,
            )
            for i in range(0, len(mzs))
        )
        record += ')'
        with open(mass_tree_file, 'w') as mass_tree:
            mass_tree.write(record)


def main():
    mgf_data_file, mass_tree_file = input().split()
    create_mass_tree(mgf_data_file, mass_tree_file)


if __name__ == '__main__':
    main()  # ../files/Magma_testing/Challenge-082.mgf ../files/Magma_testing/ms_tree
