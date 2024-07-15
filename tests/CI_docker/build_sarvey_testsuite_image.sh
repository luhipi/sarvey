#!/usr/bin/env bash
set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

context_dir="./context"
dockerfile="sarvey_ci.docker"
python_script='
version = {}
with open("../../sarvey/version.py") as version_file:
    exec(version_file.read(), version)
print(version["__version__"])
'
version=`python -c "$python_script"`
tag="sarvey_ci:$version"
gitlab_runner="sarvey_gitlab_CI_runner"

echo "#### Build runner docker image"
if [[ "$(docker images ${tag} | grep ${tag} 2> /dev/null)" != "" ]]; then
  docker rmi ${tag}
fi
DOCKER_BUILDKIT=1 docker build ${context_dir} \
    --no-cache \
    -f ${context_dir}/${dockerfile} \
    -m 20G \
    -t ${tag}
ls
