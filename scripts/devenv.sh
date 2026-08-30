# Sourced, not executed: exports the host identity that devcontainer.json
# picks up via ${localEnv:...}.
#   . scripts/devenv.sh

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    echo "devenv.sh must be sourced: . ${BASH_SOURCE[0]}" >&2
    exit 1
fi

devenv_init() {
    local name email
    name=$(git config user.name)   || return 1
    email=$(git config user.email) || return 1

    if [ -z "$name" ] || [ -z "$email" ]; then
        echo "devenv: git user.name/user.email unset on the host" >&2
        return 1
    fi

    # HOST_UID/HOST_GID rather than UID/GID: bash marks UID readonly, so
    # any later reassignment of it fails.
    export HOST_UID="$(id -u)"
    export HOST_GID="$(id -g)"
    export GIT_AUTHOR_NAME="$name"
    export GIT_AUTHOR_EMAIL="$email"
    export GIT_COMMITTER_NAME="$name"
    export GIT_COMMITTER_EMAIL="$email"
}

devenv_init || echo "devenv: environment not set" >&2
unset -f devenv_init
