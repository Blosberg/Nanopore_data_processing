export GUIX_PROFILE=${PWD}/.guix-profile
export GUIX_LOCPATH=${GUIX_PROFILE}/lib/locale


export SSL_CERT_DIR="${GUIX_PROFILE}/etc/ssl/certs"
export SSL_CERT_FILE="${GUIX_PROFILE}/etc/ssl/certs/ca-certificates.crt"
export GIT_SSL_CAINFO="${SSL_CERT_FILE}"

source ${GUIX_PROFILE}/etc/profile   # now the environmental variables will be augmented (from the front)


