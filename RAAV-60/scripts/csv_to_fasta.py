import csv

with open('LUTdna.csv', 'r') as csv_file:
    reader = csv.DictReader(csv_file)
    with open('LUTdna.fna', 'w') as fasta_file:
        for row in reader:
            fasta_file.write('>{}\n{}\n'.format(row['LUTnr'], row['Seq']))