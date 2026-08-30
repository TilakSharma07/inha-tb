#!/usr/bin/env bash
# One command to publish the InhA / M. tuberculosis pipeline to GitHub.
#
#   ./publish_to_github.sh [repo-name] [--private]
#
# Handles git identity and gh login if they are not set up yet. Stops on
# anything unexpected rather than guessing. Never force-pushes, never deletes.
set -euo pipefail

REPO="${1:-inha-tb}"
VIS="--public"
[[ "${2:-}" == "--private" ]] && VIS="--private"

cd "$(dirname "$(readlink -f "$0")")"

# 1. right folder?
for f in README.md run_all.sh src docs; do
  [[ -e "$f" ]] || { echo "ERROR: $f missing - run this from inside the inha-tb/ folder."; exit 1; }
done

# 2. git present?
command -v git >/dev/null || { echo "ERROR: git not installed.  sudo apt install git"; exit 1; }

# 3. identity - ask instead of failing, so this really is one command
if ! git config --get user.name >/dev/null 2>&1; then
  read -r -p "Your name for git commits: " GN
  git config --global user.name "$GN"
fi
if ! git config --get user.email >/dev/null 2>&1; then
  echo "Tip: use USERNAME@users.noreply.github.com to keep your real address off public commits."
  read -r -p "Your email for git commits: " GE
  git config --global user.email "$GE"
fi

# 4. init only if needed - never clobber an existing repo
[[ -d .git ]] && echo "note: existing git repo reused" || git init -b main
git add .

echo
echo "=== will commit $(git diff --cached --name-only | wc -l) files, $(du -sh --exclude=.git . | cut -f1) total ==="
git diff --cached --name-only | sed 's/^/  /' | head -25
[[ $(git diff --cached --name-only | wc -l) -gt 25 ]] && echo "  ...and more"
echo
read -r -p "Publish to github.com as '$REPO' ($VIS)? [y/N] " ok
[[ "$ok" == [yY] ]] || { echo "aborted - nothing pushed."; exit 0; }

git diff --cached --quiet || git commit -q -m "InhA / M. tuberculosis virtual screening pipeline"

# 5. push
if command -v gh >/dev/null 2>&1; then
  gh auth status >/dev/null 2>&1 || gh auth login
  gh repo create "$REPO" $VIS --source=. --remote=origin --push
  echo; echo "Done: $(gh repo view "$REPO" --json url -q .url 2>/dev/null || echo "pushed")"
else
  echo
  echo "gh CLI not installed. Either:"
  echo "  sudo apt install gh     # then re-run this script"
  echo "or create an EMPTY repo '$REPO' at https://github.com/new (no README/licence) and run:"
  echo "  git remote add origin https://github.com/USERNAME/$REPO.git"
  echo "  git push -u origin main"
fi
