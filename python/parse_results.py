"""
Quick and dirty parse function for the summary printed by uppaal.
"""
def parse_results(filename):
    data = {}
    with open(filename, 'r') as reader:
        result = reader.read()
        lines = result.split('\n')
        output = lines[-7:-2]
        for action in output:
            line = action.split(" ")
            key = f'{" ".join(line[2:-2])} [{line[-1]}]'
            data[key] = line[-2]
    return data