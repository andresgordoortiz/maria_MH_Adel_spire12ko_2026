#!/usr/bin/env bash
# Install the pinned official Linux rig binary without administrator privileges.
set -euo pipefail
prefix=${RIG_INSTALL_PREFIX:-"$HOME/.local"}
if [[ -x "$prefix/bin/rig" ]]; then
  "$prefix/bin/rig" --version
  exit 0
fi
[[ $(uname -s) == Linux ]] || { echo 'Automatic rig installation supports Linux only.' >&2; exit 1; }
case $(uname -m) in
  x86_64)
    arch=x86_64
    digest=a3ac2dd9c675247c8d5de1c8d650090df9c99e0e28e449e0546e3a54c80107cd
    ;;
  aarch64|arm64)
    arch=aarch64
    digest=46cd85e5dcbe3748c0a13e19c67333ccc9c52b07e0642a07c00cb8e7d9e6824c
    ;;
  *) echo 'Unsupported Linux architecture for rig.' >&2; exit 1 ;;
esac
stage=$(mktemp -d)
trap 'rm -rf "$stage"' EXIT
url="https://github.com/r-lib/rig/releases/download/v0.10.0/rig-linux-${arch}-0.10.0.tar.gz"
echo "Installing rig 0.10.0 into $prefix/bin"
curl --fail --location --retry 3 --connect-timeout 30 "$url" -o "$stage/rig.tar.gz"
printf '%s  %s\n' "$digest" "$stage/rig.tar.gz" | sha256sum --check --status
tar -xzf "$stage/rig.tar.gz" -C "$stage" bin/rig
"$stage/bin/rig" --version
mkdir -p "$prefix/bin"
install -m 755 "$stage/bin/rig" "$prefix/bin/rig"
