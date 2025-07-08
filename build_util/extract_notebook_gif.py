import json

from nbconvert.preprocessors import Preprocessor
import os
import base64

flag = False

class ExtractGifPreprocessor(Preprocessor):
    """A preprocessor which extracts gifs.

    Because nbconvert does not do it.

    (Thanks ChatGPT!)

    """

    def preprocess_cell(self, cell, resources, cell_index):
        if 'outputs' not in cell:
            return cell, resources

        # for debugging.
        global flag
        if not flag:
            for key, val in resources.items():
                print(f"{key}: {val}")
            print()
            flag = True

        notebook_name = resources["unique_key"]  # can't fly without it
        notebook_dir = os.path.dirname(notebook_name)

        output_dir = resources.get('output_files_dir', 'notebook_files')
        output_dir_relative = os.path.relpath(output_dir, notebook_dir)

        print(notebook_name, notebook_dir, output_dir, output_dir_relative)

        os.makedirs(output_dir, exist_ok=True)

        for i, output in enumerate(cell['outputs']):
            data = output.get('data', {})
            gif_data = data.get('image/gif')
            if gif_data:
                filename = f"output_{cell_index}_{i}.gif"
                filepath = os.path.join(output_dir, filename)
                with open(filepath, 'wb') as f:
                    f.write(base64.b64decode(gif_data))
                # Replace output with HTML <img> tag
                output['data'] = {
                    'text/html': f'<img src="{output_dir_relative}/{filename}">'
                }
        return cell, resources
