""" A simple script to copy docker images specified in a json file to CPG's AR"""

import subprocess
import json
from urllib.request import urlopen
import click


@click.command()
@click.option('--docker_json_url')
def main(docker_json_url):
    """ Copies each docker image """

    if not docker_json_url:
        github_url = 'https://raw.githubusercontent.com/'
        path_to_docker_json = (
            'populationgenomics/gatk-sv/main/input_values/dockers.json'
        )
        docker_json_url = github_url + path_to_docker_json

    subprocess.run(
        ['gcloud', 'auth', 'configure-docker', 'australia-southeast1-docker.pkg.dev'],
        check=True,
    )

    response = urlopen(docker_json_url)
    docker_json = json.loads(response.read())
    docker_json.pop('name')

    cpg_docker_dict = {}
    for key in docker_json:
        initial_path = 'docker://' + docker_json[key]  # noqa: F841
        docker_image = docker_json[key].split('/')[-1]
        cpg_ar_path = (
            'australia-southeast1-docker.pkg.dev/cpg-common/images/sv/' + docker_image
        )
        destination_path = 'docker://' + cpg_ar_path  # noqa: F841
        cpg_docker_dict[key] = cpg_ar_path
        # subprocess.run(['skopeo', 'copy', initial_path, destination_path], check='True')

    with open('pipelines/gatk_sv/dockers.json', 'w') as output_json:
        json.dump(cpg_docker_dict, output_json, indent=4)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
