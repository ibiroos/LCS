import json
from collections import OrderedDict


def read_input(input_file, input_keys):
    try:
        with open(input_file, 'r') as f:
            return json.load(f, object_pairs_hook=OrderedDict)
    except FileNotFoundError:
        print(f'File not found: {input_file} ')
        if input('Do you want to create one (y/n)?') == 'n':
            quit()

        print(f'A {input_file} will be created with the next keys:\n')
        json_obj = {}
        for input_key in input_keys:
            json_obj[input_key] = input(f'key: {input_key}?\n')
        print(f'Writing a json_file: {input_file} with the next content:')
        print(json.dumps(json_obj, indent=4))
        with open(input_file, 'w') as json_file:
            json.dump(json_obj, json_file, indent=4)
        print('Done!\n')
        return read_input(input_file, input_keys)
