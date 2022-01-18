import os
import shutil


def parse_from_mgf(s):
    res = []
    cur = ''
    in_brackets = False
    for a in s:
        if (a == ',') and (not in_brackets):
            res.append(cur)
            cur = ''
        elif a == '\"':
            in_brackets = not in_brackets
        else:
            cur += a
    if cur != '':
        res.append(cur)
    return res


def _transform_from_mgf_to_mass_tree(
        from_file,
        to_file,
):
    with open(from_file) as mgf:
        mzs = []
        norm_intensities = []
        mass = None
        for line in mgf.readlines():
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
        with open(to_file, 'w') as mass_tree:
            mass_tree.write(record)


def _transform_from_csv_to_sql(
        from_file,
        to_file,
):
    pass


def transform_from_to(
        from_format,
        to_format,
        from_file,
        to_file,
):
    if from_format == to_format:
        shutil.copyfile(from_file, to_file)
    elif (from_format == 'mgf') and (to_format == 'tree'):
        _transform_from_mgf_to_mass_tree(from_file, to_file)
    elif (from_format == 'csv') and (to_format == 'db'):
        _transform_from_csv_to_sql(from_file, to_file)
    else:
        raise NotImplementedError


def get_name(conn, user_id):
    c = conn.cursor()
    c.execute("SELECT name FROM userList WHERE userID = ?", (user_id,))
    result = c.fetchone()
    if result:
        return result[0]
