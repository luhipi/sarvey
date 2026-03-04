#!/bin/bash
set -e

# Check for necessary environment variables
if [ -z "$RUNNER_TOKEN" ] || [ -z "$RUNNER_REPO_URL" ]; then
  echo "Error: RUNNER_TOKEN and RUNNER_REPO_URL environment variables must be set."
  exit 1
fi

# Configure the runner if it hasn’t been configured already
if [ ! -f ".runner" ]; then
  echo "Configuring the GitHub Actions runner..."
  ./config.sh --unattended \
    --replace \
    --labels self-hosted \
    --url "$RUNNER_REPO_URL" \
    --token "$RUNNER_TOKEN" \
    --name "${RUNNER_NAME:-docker-runner}" || exit 1
fi

# Trap SIGTERM and deregister the runner on shutdown if needed
trap './config.sh remove --token "$RUNNER_TOKEN"; exit 0' SIGTERM

echo "Starting the GitHub Actions runner..."
exec ./run.sh

