#!/usr/bin/env bash
# Render docs/overview.html (editable source) to docs/overview.png (embedded in README).
set -euo pipefail
cd "$(dirname "$0")/.."
CHROME="${CHROME:-/Applications/Google Chrome.app/Contents/MacOS/Google Chrome}"
command -v google-chrome >/dev/null && CHROME=google-chrome
"$CHROME" --headless=new --disable-gpu --hide-scrollbars --force-device-scale-factor=2 \
  --window-size=1400,1180 --screenshot="$PWD/docs/overview.png" "file://$PWD/docs/overview.html"
echo "wrote docs/overview.png"
