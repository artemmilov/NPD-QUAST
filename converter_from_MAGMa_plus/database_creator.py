import sqlite3
import zlib

from magma_plus_mol import MagmaPlusMol, MagmaPlusInitException


def create_database(csv_database_file, converted_database_file):
    conn = sqlite3.connect(converted_database_file)
    c = conn.cursor()
    try:
        c.execute(
            """CREATE TABLE molecules (id TEXT PRIMARY KEY,
                                             mim INTEGER NOT NULL,
                                             charge INTEGER NOT NULL,
                                             natoms INTEGER NOT NULL,
                                             molblock TEXT,
                                             inchikey TEXT,
                                             molform TEXT,
                                             name TEXT,
                                             reference TEXT,
                                             logp INT)""",
        )
        conn.commit()
        print("{0}.db created".format(converted_database_file))
    except:
        print(
            "{0} already exists (or error creating it)".format(
                converted_database_file,
            ),
        )
        return

    with open(csv_database_file) as csv_database:
        ok = 0
        not_ok = 0
        for line in csv_database.readlines()[1:]:
            if line != '':
                try:
                    magma_plus_mol = MagmaPlusMol(line)
                    c.execute(
                        '''INSERT INTO molecules (id, mim, charge, natoms, molblock, inchikey,
                         molform, name, reference, logp) VALUES (?,?,?,?,?,?,?,?,?,?)''',
                        (
                            magma_plus_mol.scan,
                            int(round(magma_plus_mol.mass * 1e6)),
                            magma_plus_mol.charge,
                            magma_plus_mol.n_atoms,
                            sqlite3.Binary(zlib.compress(magma_plus_mol.mol_block.encode('utf-8'))),
                            magma_plus_mol.inchi_key,
                            magma_plus_mol.mol_form,
                            magma_plus_mol.name,
                            magma_plus_mol.reference,
                            int(round(magma_plus_mol.logp * 10)),
                        ),
                    )
                    ok += 1
                except MagmaPlusInitException:
                    not_ok += 1
    conn.commit()

    print("Creating index ...")
    c.execute('PRAGMA temp_store = 2')
    c.execute(
        'CREATE INDEX idx_cover ON molecules (charge, mim, natoms, reference, molform, inchikey, name, molblock, logp)',
    )
    conn.commit()
    print('Ok: {0}, not ok: {1}.'.format(ok, not_ok))


def main():
    csv_database_file, converted_database_file = input().split()
    create_database(csv_database_file, converted_database_file)


if __name__ == '__main__':
    main()  # ../files/Magma_testing/Challenge-082.csv ../files/Magma_testing/data.db
    # python MAGMa_plus.py annotate -c 0 -d 0 -b 3 -w 1 ->
    # --db_options ../NPD-QUAST/Magma_testing/data.db ../NPD-QUAST/Magma_testing/output.db
