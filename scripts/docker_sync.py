# Syncing Docker Images

import subprocess

from urllib.request import urlopen

import json


def main():
    """ Move each docker image specified in json """

    subprocess.run(
        ["gcloud", "auth", "configure-docker", "australia-southeast1-docker.pkg.dev"]
    )
    # Eventually  create an optional command line input to parse any json

    docker_json_url = "https://raw.githubusercontent.com/populationgenomics/gatk-sv/main/input_values/dockers.json"

    response = urlopen(docker_json_url)
    docker_json = json.loads(response.read())
    docker_json.pop("name")

    for key in docker_json:
        initial_path = "docker://" + docker_json[key]
        docker_image = docker_json[key].split("/")[-1]
        destination_path = (
            "docker://australia-southeast1-docker.pkg.dev/cpg-common/images/sv/"
            + docker_image
        )
        subprocess.run(["skopeo", "copy", initial_path, destination_path])


if __name__ == "__main__":
    main()