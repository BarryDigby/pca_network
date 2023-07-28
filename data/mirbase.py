def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        prev_key = None
        prev_values = []

        for line in infile:
            key, values = line.strip().split('\t')
            values = values.split(';')

            if prev_key is None:
                prev_key = key
                prev_values = values
            elif key != prev_key:
                write_output(outfile, prev_key, prev_values)
                prev_key = key
                prev_values = values
            else:
                prev_values.extend(values)

        write_output(outfile, prev_key, prev_values)

def write_output(outfile, key, values):
    for value in values:
        output_line = f"{key}\t{value}\n"
        outfile.write(output_line)

if __name__ == "__main__":
    input_file_path = "mirbase_aliases.txt"
    output_file_path = "output.txt"
    process_file(input_file_path, output_file_path)

