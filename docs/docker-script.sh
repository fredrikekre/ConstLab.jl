#!/bin/bash

JULIA_VERSION=$1
JULIABIN=/test/julia-${JULIA_VERSION}/bin/julia

echo "Installing curl"
apt-get update -qq && apt-get install --no-install-recommends -qq curl

echo "Installing Julia"
/test/install-julia.sh $JULIA_VERSION

echo "Installing documentation dependencies"
$JULIABIN --color=yes --project=docs -e "import Pkg; Pkg.instantiate();"

echo "Building and deploying documentation"
$JULIABIN --color=yes --project=docs docs/make.jl
