from flask import Flask, jsonify, request, abort
from flask_cors import CORS

import json
from os import remove
from tempfile import NamedTemporaryFile
import subprocess

application = Flask(__name__)
CORS(application)


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


@application.route('/api/singlepoint', methods=['POST'])
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

@application.route('/api/matlab', methods=['POST'])
def matlab_test():
    if not request.json:
        abort(400)

    with NamedTemporaryFile(delete=False, suffix='.json', mode='wt') as inf:
        with NamedTemporaryFile(delete=False, suffix='.json') as outf:
            in_fn = inf.name
            out_fn = outf.name
            # write request json to temporary file for matlab to read
            inf.write(json.dumps(request.json))

    # Create command line args, or incorporate the SESAM function?
    # sesam_command = '/home/ec2-user/code_kim/VECTOR/api/vector_calc.sh'
    sesam_command = '/opt/python/current/app/vector_calc.sh'

    # Call matlab SESAM
    result = subprocess.run([sesam_command, in_fn, out_fn], check=True)

    # TODO handle errors
    #print(result.returncode)

    with open(out_fn, "r") as f:
        output_payload = f.read()

    remove(in_fn)
    remove(out_fn)

    # print(output_payload)

    output_dict = json.loads(output_payload)
    return jsonify(output_dict)


if __name__ == '__main__':
    application.run(debug=True)
