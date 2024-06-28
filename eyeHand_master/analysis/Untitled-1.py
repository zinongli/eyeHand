def extract_string_rows(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.strip() and not line.strip().split()[0].isdigit():
                outfile.write(line)

# Replace 'input.txt' with the name of your input file and 'output.txt' with the desired output file name
input_file = 'pilo0425.txt'
output_file = '0425Event.txt'
extract_string_rows(input_file, output_file)

def parse_trials(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        trial_number = None
        saccade_onset_found = False
        for line in infile:
            parts = line.strip().split()
            if parts[0] == 'MSG' and 'TRIAL_START' in parts:
                trial_number = parts[-1]
                saccade_onset_found = False
            elif parts[0] == 'MSG' and 'SaccadeOnset' in parts:
                saccade_onset_found = True
            elif parts[0] == 'ESACC' and saccade_onset_found:
                outfile.write(f"{trial_number},{','.join(parts[4:])}\n")
                saccade_onset_found = False

# Replace 'input.txt' with the name of your input file and 'output.txt' with the desired output file name
input_file = '0425Event.txt'
output_file = '0425Sacc1st.csv'
parse_trials(input_file, output_file)
