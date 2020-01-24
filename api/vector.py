from flask import Flask, jsonify, request, abort
from flask_cors import CORS

import subprocess

app = Flask(__name__)
CORS(app)


def process_input_payload(data):
    """
    Process the payload from the frontend and properly format
    a the sesam code command line call.
    """
    # TODO: Read the payload and update the command
    print(data)

    # cmd = 'matlab' + command line args from parsing json data

    # Return a simple function call for now
    return 'ls -la'


def process_sesam_output(data):
    """Process the SESAM code output to properly format a return payload."""
    # TODO: Update the payload based on what is in data
    payload = {"dragCoefficient": 2.8,
               "projectedArea": 3,
               "forceCoefficient": 1.2}
    return payload


@app.route('/api/singlepoint', methods=['POST'])
def sesam_run():
    if not request.json:
        abort(400)

    # Create command line args, or incorporate the SESAM function?
    sesam_command = process_input_payload(request.json)

    # Call matlab SESAM
    output = subprocess.Popen(sesam_command, shell=True,
                              stdout=subprocess.PIPE).stdout.read()

    print(output)

    # Process SESAM output to properly format it on return
    output_payload = process_sesam_output(output)

    return jsonify(output_payload)


if __name__ == '__main__':
    app.run(debug=True)
