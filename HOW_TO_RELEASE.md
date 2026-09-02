# How to release

The document describes how to release `gpsea` to *PyPi*.

## Release checklist

- update documentation and notebooks to present the latest features 
- create and checkout a release branch
- remove deprecated methods targeted for removal in this version. The `TODO` markers are labeled using 
  the target version (e.g. `TODO[v0.3.0]`)
- bump versions to a release:
  - `pyproject.toml`
- ensure the CI passes
- merge to `main`
- create a GitHub release from the latest `main` commit, including a new tag. The `release.yml` workflow takes over and deploys the new release to PyPi.
- merge `main` to `develop`
- bump versions to a `dev` version to begin next development iteration

